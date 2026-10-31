"""Grid, checkpoint and renderer checks independent of the long experiment."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from copa2026 import figures


def full_rows():
    rows = []
    for dataset in figures.DATASET_NAMES:
        for k in range(1, 101):
            for value in figures.DEFAULT_LAMBDAS:
                rows.append(dict(Dataset=dataset, Implementation="normalized", K=k, Lambda=value,
                                 NCM="ellipsoid_pdf", Epsilon=100, DetectionRate=k/100,
                                 FalseAlarmRate=.01))
            rows.append(dict(Dataset=dataset, Implementation="normalized", K=k, Lambda=.001,
                             NCM="knn_distance", Epsilon=100, DetectionRate=.8,
                             FalseAlarmRate=.02))
    return rows


def write_rows(path, rows):
    with path.open("w", newline="", encoding="utf8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_default_plan_is_one_joint_grid():
    jobs=figures._jobs(figures.DATASET_NAMES,range(1,101),figures.DEFAULT_LAMBDAS,True,100.)
    assert len(jobs)==2400
    assert [job["Number"] for job in jobs]==list(range(1,2401))
    assert sum(job["NCM"]=="knn_distance" for job in jobs)==400
    assert sum(job["NCM"]=="ellipsoid_pdf" and job["Lambda"]==.001 for job in jobs)==400


@pytest.mark.parametrize("fault",["missing","duplicate","epsilon","lambda","rate","k"])
def test_renderer_rejects_invalid_full_grid(tmp_path,fault):
    rows=full_rows()
    if fault=="missing": rows.pop()
    elif fault=="duplicate": rows[-1]=rows[0].copy()
    elif fault=="epsilon": rows[0]["Epsilon"]=10
    elif fault=="lambda": rows[0]["Lambda"]=1.001e-5
    elif fault=="rate": rows[0]["DetectionRate"]=float("nan")
    else: rows[0]["K"]=1.5
    file=tmp_path/"Performance.csv"
    write_rows(file,rows)
    with pytest.raises(ValueError): figures._validated_plot_rows(file)


def test_complete_renderer_keeps_all_curves_and_outputs(tmp_path,monkeypatch):
    root=tmp_path/"repository"/"python"
    (root/"results").mkdir(parents=True)
    monkeypatch.setattr(figures,"PACKAGE_ROOT",root)
    source=root/"results"/"Performance.csv"
    write_rows(source,full_rows())
    output=figures.reconstruct_copa_figures(source,export_pdf=True)
    metadata=json.loads((output/"reconstruction_metadata.json").read_text())
    assert metadata["RowsUsed"]==2400
    assert metadata["ExperimentsRecomputedByRenderer"] is False
    assert [record["CurveCount"] for record in metadata["PlotRecords"]]==[2]*8+[5]*8
    assert all(record["PointsPerCurve"]==100 for record in metadata["PlotRecords"])
    assert len(list(csv.DictReader((output/"PerformanceUsed.csv").open())))==2400
    for name in metadata["ArtifactNames"]:
        assert (output/name).stat().st_size>100
        assert "normalized" not in name


def test_resume_checks_sources_inputs_and_prepared_bytes(tmp_path,monkeypatch):
    root=tmp_path/"repository"/"python"
    data=root/"datasets"
    data.mkdir(parents=True)
    for name,file in figures.FILES.items():
        np.savez_compressed(data/file,X=np.arange(40,dtype=float).reshape(20,2))
    monkeypatch.setattr(figures,"PACKAGE_ROOT",root)
    monkeypatch.setattr(figures,"_source_hashes",lambda:{"fixture":"unchanged"})
    options=dict(data_dir=data,dataset_names=["protein"],k_list=[1],lambda_list=[.001],
                 reps=1,max_jobs=0,workers=1)
    rows,directory=figures.run_full_figure_sweeps(**options)
    assert rows==[]
    assert json.loads((directory/"progress.json").read_text())["Status"]=="partial"
    assert len(list((directory/"prepared").glob("*.npy")))==4
    _,same=figures.run_full_figure_sweeps(**dict(options,workers=16,output_root=directory))
    assert same==directory
    with pytest.raises(ValueError,match="Resume rejected"):
        figures.run_full_figure_sweeps(**dict(options,reps=2,output_root=directory))
    scipy_version=figures.scipy.__version__
    monkeypatch.setattr(figures.scipy,"__version__","changed-runtime")
    with pytest.raises(ValueError,match="Resume rejected"):
        figures.run_full_figure_sweeps(**dict(options,output_root=directory))
    monkeypatch.setattr(figures.scipy,"__version__",scipy_version)
    monkeypatch.setattr(figures,"_source_hashes",lambda:{"fixture":"changed"})
    with pytest.raises(ValueError,match="Resume rejected"):
        figures.run_full_figure_sweeps(**dict(options,output_root=directory))
    monkeypatch.setattr(figures,"_source_hashes",lambda:{"fixture":"unchanged"})
    np.save(directory/"prepared"/"protein.npy",np.zeros((20,2)))
    with pytest.raises(ValueError,match="prepared data changed"):
        figures.run_full_figure_sweeps(**dict(options,output_root=directory))


def test_explicit_output_parent_is_supported(tmp_path,monkeypatch):
    monkeypatch.setattr(figures,"PACKAGE_ROOT",tmp_path/"repository"/"python")
    output=figures.new_output_dir("run",tmp_path/"custom_outputs")
    assert output.is_dir()
    assert output.parent==tmp_path/"custom_outputs"

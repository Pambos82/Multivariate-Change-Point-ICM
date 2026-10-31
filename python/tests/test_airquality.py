"""Air Quality preparation and stream behavior using synthetic local data."""
from datetime import datetime,timedelta
import json

import numpy as np
import pytest
from openpyxl import Workbook

from copa2026 import figures


def make_spreadsheet(path):
    book=Workbook(write_only=True)
    sheet=book.create_sheet()
    sheet.append(["Date","Time","CO(GT)","PT08.S1(CO)","NMHC(GT)","C6H6(GT)","PT08.S2(NMHC)",
                  "NOx(GT)","PT08.S3(NOx)","NO2(GT)","PT08.S4(NO2)","PT08.S5(O3)","T","RH","AH"])
    for row in range(9357):
        measurements=[float(i)+.125 for i in range(13)]
        measurements[2]=-200  # Dropped NMHC must not remove otherwise complete rows.
        if row>=6941: measurements[0]=-200
        sheet.append([datetime(2004,3,10)+timedelta(days=row//24),(row%24)/24,*measurements])
    book.save(path)


def test_original_xlsx_schema_precision_cleaning_and_week_index(tmp_path,monkeypatch):
    root=tmp_path/"repository"/"python"
    root.mkdir(parents=True)
    monkeypatch.setattr(figures,"PACKAGE_ROOT",root)
    source=tmp_path/"AirQualityUCI.xlsx"
    make_spreadsheet(source)
    destination=figures.prepare_airquality_input(source,tmp_path/"local_air")
    with np.load(destination,allow_pickle=False) as stored:
        assert stored["X"].shape==(6941,12)
        assert stored["X"][0,0]==.125
        assert stored["X"][0,2]==3.125
        assert stored["week"][0]==1
        assert stored["week"][168]==2
        assert stored["dates"][0]==stored["dates"][23]
    with pytest.raises(FileExistsError): figures.prepare_airquality_input(source,destination.parent)
    inside_fixture=figures.prepare_airquality_input(source,root/"results"/"synthetic_air_input")
    assert inside_fixture.is_file()
    with pytest.raises(ValueError,match="XLSX|xlsx"):
        figures.prepare_airquality_input(tmp_path/"rounded.csv",tmp_path/"csv_input")


def test_stream_uses_log_alarm_and_shared_uniforms(tmp_path):
    x=np.random.RandomState(3).normal(size=(32,3))
    week=np.arange(32)//7+1
    uniforms=np.linspace(.03,.97,12)
    stream=figures.run_real_stream(x,week,TR_size=20,K=5,uniforms=uniforms)
    hits=np.flatnonzero(stream["logM"]>np.log(100))
    assert stream["alarm_index"]==(int(hits[0]+1) if len(hits) else None)
    np.testing.assert_array_equal(stream["week_test"],week[20:])
    repeat=figures.run_real_stream(x,week,TR_size=20,K=5,uniforms=uniforms)
    np.testing.assert_array_equal(stream["pval"],repeat["pval"])
    np.testing.assert_array_equal(stream["diagnostics"]["selected"],repeat["diagnostics"]["selected"])


def test_run_real_writes_complete_figure3_outputs(tmp_path,monkeypatch):
    root=tmp_path/"repository"/"python"
    root.mkdir(parents=True)
    monkeypatch.setattr(figures,"PACKAGE_ROOT",root)
    air=tmp_path/"local_air"
    air.mkdir()
    x=np.random.RandomState(4).normal(size=(3010,12))
    dates=np.datetime64("2004-03-10")+np.arange(len(x)).astype("timedelta64[D]")
    week=np.arange(len(x))//7+1
    np.savez_compressed(air/"airquality.npz",X=x,week=week,dates=dates.astype("U10"))
    stream,summary,output=figures.run_real(data_dir=air)
    assert summary["MonitoredSamples"]==10
    assert summary["FiniteLogM"] is True
    assert len(stream["rowTimes"])==10
    for name in ["figure3.png","figure3.pdf","airquality_results.npz","airquality_results.md",
                 "Summary.csv","SelectionFrequencies.csv","NormalizationConstants.csv","provenance.json"]:
        assert (output/name).is_file()
    assert "normalized" not in (output/"airquality_results.md").read_text().lower()
    provenance=json.loads((output/"provenance.json").read_text())
    assert provenance["ReferenceSize"]==3000
    assert provenance["NumPyVersion"]==np.__version__
    assert provenance["SciPyVersion"]==figures.scipy.__version__

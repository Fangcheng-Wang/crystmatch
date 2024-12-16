import os
import glob
import shutil
import subprocess

sep = os.path.sep

def clean_ext(ext_list):
    for ext in ext_list:
        for file in glob.glob("*." + ext):
            os.remove(file)

def clean_csm():
    for item in os.listdir("."):
        if os.path.isdir(item) and item.startswith("CSM_"):
            shutil.rmtree(item)

def test_version():
    result = subprocess.run(
        ["crystmatch", "-v"], 
        capture_output=True, text=True, timeout=10
    )
    assert result.returncode == 0
    assert "crystmatch" in result.stdout

def test_help():
    result = subprocess.run(
        ["crystmatch", "-h"], 
        capture_output=True, text=True, timeout=5
    )
    assert result.returncode == 0
    assert "wfc@pku.edu.cn" in result.stdout

def test_enumeration():
    clean_ext(["npz", "pdf", "csv"])
    result = subprocess.run(
        ["crystmatch"], input=f"enumeration\n3\n0.1\ntest{sep}wz.txt\ntest{sep}zb.txt\n",
        capture_output=True, text=True, timeout=10
    )
    assert result.returncode == 0
    assert "Found 0" in result.stdout
    assert "Found 2" in result.stdout
    assert "Found 10" in result.stdout
    result = subprocess.run(
        ["crystmatch", "-E", "3", "0.1", "-I", f"test{sep}wz.txt", "-F", f"test{sep}zb.txt"],
        capture_output=True, text=True, timeout=10
    )
    assert result.returncode == 0
    assert "A total of 12" in result.stdout
    assert "'CSM_LIST-m3s0.10-1.npz'" in result.stdout
    assert os.path.exists("CSM_LIST-m3s0.10-1.npz")
    assert os.path.exists("PLOT-m3s0.10-1.pdf")
    assert os.path.exists("TABLE-m3s0.10-1.csv")

def test_read():
    clean_csm()
    result = subprocess.run(
        ["crystmatch"], input=f"read\nCSM_LIST-m3s0.10.npz\n7 11\n",
        capture_output=True, text=True, timeout=5
    )
    assert result.returncode == 0
    result = subprocess.run(
        ["crystmatch", "-R", "CSM_LIST-m3s0.10.npz", "-e", "7", "11"],
        capture_output=True, text=True, timeout=5
    )
    assert result.returncode == 0
    assert os.path.exists(f"CSM_11-1{sep}POSCAR_I")

def test_single():
    result = subprocess.run(
        ["crystmatch"], input=f"single\nCSM_7{sep}POSCAR_I\nCSM_7{sep}POSCAR_F\n",
        capture_output=True, text=True, timeout=5
    )
    assert result.returncode == 0
    assert "rmss = 0.0842" in result.stdout
    assert "rmsd = 1.1126" in result.stdout
    result = subprocess.run(
        ["crystmatch", "-S", "-I", f"CSM_11{sep}POSCAR_I", "-F", f"CSM_11{sep}POSCAR_F"],
        capture_output=True, text=True, timeout=5
    )
    assert result.returncode == 0
    assert os.path.exists(f"CSM_single-1{sep}POSCAR_I")
    assert os.path.exists(f"CSM_single-1{sep}POSCAR_F-optimized")
    
def test_clean():
    clean_ext(["npz", "pdf", "csv"])
    clean_csm()
    assert 1 + 1 == 2
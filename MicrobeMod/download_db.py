import subprocess


print("Downloading database from FigShare...")

cmd = "wget https://figshare.com/ndownloader/files/42960463 -O db.tar.gz"
output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()

print("Extracting database")

cmd = "tar -xvf db.tar.gz"
output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()

output = subprocess.check_output(
    "rm db.tar.gz", shell=True, stderr=subprocess.STDOUT
).decode()
print("Database complete!")

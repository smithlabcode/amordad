Amordad
=======
# NOTE
We need to think of the way we are defining our config files. There is loads of duplication
and loads of ‘reinventing the format wheel’. This needs to be simply chucked out.
Above all, all the formats should be documented!

# Setting it all up

1. Create a config file: `config.txt`
```
/media/data1/Development_Version_Controlled/Research/wenzhenl-amordad/engine/mgrast_mgmap.txt
/media/data1/Development_Version_Controlled/Research/wenzhenl-amordad/engine/qhf.txt
/media/data1/Development_Version_Controlled/Research/wenzhenl-amordad/engine/qht.txt
/media/data1/Development_Version_Controlled/Research/wenzhenl-amordad/engine/graph/graph_mgrast348.out
```
2. Download these four files. I took a ‘shortcut’ and downloaded EVERYTHING:
```
rsync -avh skchoudh@smithdb.usc.edu:/home/wenzhenl/amordad  /media/data1/Development_Version_Controlled/Research/wenzhenl-amordad
```

3. Download MGRAST database from server

```
rsync -avh skchoudh@hpc-login3.usc.edu:/home/rcf-47/andrewds/panfs/Amordad_data/metagenomes/MGRAST mgrast
```

4. Edit `mgrast_mgmap.txt` to reflect the absolute path of k-mer count files poiting to location you got from step 3

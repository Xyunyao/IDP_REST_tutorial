# Installing gromacs patched with Plumed2

install.sh has a working script to install gromacs 2024.3 and plumed2 v2.9 in the users $HOME folder under $HOME/opt.
From this directory simply perform the two commands to install. 
```bash
chmod +x install.sh
sh install.sh
```

Feel free to modify the installation directory.
Be aware, if you change the installation directory to a location requiring access above the users privileges the script will need to be performed via sudo. 

```bash
chmod +x install.sh
sudo sh ./install.sh
```

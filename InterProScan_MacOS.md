# InterProScan install on MacOS (tested on M1,2,3)
## Set-up the lima shell
### Install lima and docker
There is no version of InterProScan for macOS, therefore it has to be run in a docker container. For this, you have to have docker installed (both in the command line, with `brew` and Docker Desktop https://docs.docker.com/desktop/setup/install/mac-install/), and lima (as a virtual linux shell, also via `brew`). 

```bash
brew install lima
brew install docker
## Test if they are installed
docker -h
limactl shell --help
```

### Modify writability in lima.yaml config
#### Short version
Find the following section in the lima config file (default location: "~/.lima/apptainer/lima.yaml"): 

```yaml
mounts:
- location: "~"

- location: "{{.GlobalTempDir}}/lima"
  mountPoint: /tmp/lima
  writable: true
```

Replace with: 

```yaml
mounts:
- location: "/Users/sinaedaa"
  writable: true

- location: "{{.GlobalTempDir}}/lima"
  mountPoint: /tmp/lima
  writable: true
```

#### Long explanation
Lima also includes `apptainer` (newer instance of Singularity), which is what Snakemake will use to run InterProScan. By default however, a lima installation on macOS sets the `$HOME` directory as unwriteable, which also affects all down-branching directories. Additionally, there are some incompatibilities between the default accessible memory and cpus given to the lima virtual machine, and what InterProScan uses by default. I changed the lima configuration file to provide more memory to the VM, otherwise the pipeline will automatically fail, as it runs out of memory. 

The lima config file should be located (by default), here `~/.lima/apptainer/lima.yaml`. Simply open it, and change some things.

All the way down the file, you will find the `mounts` section, which should look something like this:

```yaml
mounts:
- location: "~"

- location: "{{.GlobalTempDir}}/lima"
  mountPoint: /tmp/lima
  writable: true
```

By default, when you open a lima shell (`limactl shell apptainer`), it will open it in the directory you were in when you opened it. If this directory is anywhere inside your $HOME, then you will not be able to write anything there. Your $HOME should correspond to `/Users/yourcomputername`, and on the host computer, it also corresponds to "~". The host `~` and the lima `~` are not the same, and the `~` that is specified in the lima.yaml file, refers to the lima `~`. I changed the first mount location to my host $HOME (`/Users/myname`), and added "writable" below, so that lima gets access to your $HOME directory and can work inside the snakemake pipeline. 

```yaml
mounts:
- location: "/Users/sinaedaa"
  writable: true

- location: "{{.GlobalTempDir}}/lima"
  mountPoint: /tmp/lima
  writable: true
```
### Tell lima to use Rosetta
In order to correctly simulate a Linux operating system with `apptainer`, we need to use Rosetta on macOS. You don't need to install anything, just add these lines to `lima.yaml`:

```yaml
## For Rosetta
minimumLimaVersion: 1.1.0
vmType: "vz"
rosetta:
  # Enable Rosetta for Linux.
  # Hint: try `softwareupdate --install-rosetta` if Lima gets stuck at `Installing rosetta...`
  enabled: true
  # Register rosetta to /proc/sys/fs/binfmt_misc
  binfmt: true
```

### Adapt VM memory usage
Finally, in order to provide more working memory to the lima shell, add this to `lima.yaml`:

```yaml
## For CPU usage
cpus: 4
memory: "16GiB"
disk: "100GiB"
```

### Restart lima to apply changes
If lima was running in the background (likely), after changing the lima.yaml file, run this to stop and restart the shell with the new parameters: 

```bash
limactl stop apptainer
limactl start apptainer
```

### Create a symlink to apptainer.lima
By default, the apptainer instance installed with lima is called `apptainer.lima`, which snakemake is unable to find. Snakemake looks for **singularity** and/or **apptainer**, present in your `$PATH`. 

```bash
## Look where your apptainer.lima is located
which apptainer.lima # for me: /opt/homebrew/bin/apptainer.lima
## Make symlinks to a directory that is present in your path
echo $PATH # I have one dir called: $HOME/.local/bin/
ln -s /opt/homebrew/bin/apptainer.lima $HOME/.local/bin/apptainer
ln -s /opt/homebrew/bin/apptainer.lima $HOME/.local/bin/singularity
```

I added a link to `singularity` as well, as my personal version of snakemake always first looks for Singularity, and returns an error if it cannot find it.

### Testing docker and apptainer
As a test, you can try both of these commands, which should return no error. **Note**: the second one will return a "will not overwrite" error if you run it a second time and the sif image was already created. 
```bash
docker pull interpro/interproscan:5.75-106.0
apptainer pull docker://interpro/interproscan:5.75-106.0

# remove the created interproscan sif file, we don't need it, as it will be created during the workflow, and stored here .snakemake/singularity/
rm -rf interproscan_5.75-106.0.sif 
```

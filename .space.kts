job("test-yaml-files") {
	container(displayName = "mamba install yaml", image="condaforge/mambaforge:latest") {
        mountDir  = "/mnt/space"
        workDir = "/mnt/space/work"
    	shellScript {
            interpreter = "/bin/bash"
           	location = "/mnt/space/work/scn-pipeline/automation/mamba_create.sh"
        }
    }
}
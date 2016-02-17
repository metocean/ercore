Wrapper for ERI operations
===

This is the code for running ERI domains in for MOV, to be running using the Scheduler wrapper API.

Workflow
===


 * Request data from UDS
 * Create configuration files based on sites from MOV and defined boundery
 * Run the model
 * Generate Kernel Density blobs with the particle clouds to be saved in json


Thoughts
===

 - Create ctl files likewise others models in the house
 - Remove dependency from msl_actions


Testing
===

Go into tests folder and start do a:

```
nosetests -sv
```

Fell free to add more tests as needed.






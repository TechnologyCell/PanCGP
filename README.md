# PanCGP:Pan-genome and Comparative Genome analysis Pipeline tool

PanCGP (Pan-genome and Comparative Genome analysis Pipeline) tool is designed for microbiologist having minimum knowledge of computational analysis. The objective is to facilitate the users to perform computational pan-genomic analysis in a fraction of time with just a list of sequenced proteins as compared to other software packages already available for the purpose. It is designed to automatically scale depending on the available resources ranging from a personal computer to a hybrid cluster as well as the size of input. It also addresses the issue of large amount of data processing for comparative analysis of multiple proteins.

## Instalation guide
Download and install Ubuntu, Ubuntu 12.04 or latter but Ubuntu 12.04 desktop version is recommended because all tools are tested at this version. You may also use latter version but you should take care of compatibilities of all relevant tools.
- Open terminal using ALT+CTRL+T and run following commands one by one
```sh
    $ sudo apt-get update
	$ sudo apt-get install gcc build-essential
	$ sudo apt-get install default-jre
	$ sudo apt-get install default-jdk
	$ sudo apt-get install -y autotools-dev g++ build-essential openmpi1.5-bin openmpi1.5-doc libopenmpi1.5-dev
	$ sudo apt-get autoremove
	$ sudo apt-get install perl
	$ sudo apt-get install bioperl
```
    (At any step, if system prompts for any input, please press Enter and use default option.)
- Download $mpjexpress$ library from here. (Current latest is 0.44)
- Once downloaded, extract it to a directory (suppose that you have downloaded the it and extracted in your home directory then the path will be  _/home/< username>/mpj/mpj-v\_044_)
- Now we need to set environment variables for mpj express. These variables can be added in .bashrc file. The variables need to be added in start of file otherwise may cause problem.
	- Open .bashrc: ```$ vi .bashrc```
	- For java path: ```$ export JAVA_HOME=/usr/java/latest```
	- For mpj express path: ```$ export MPJ_HOME=/home/<username>/mpj/mpj-v_044```
	- Add to path variable: ```$export PATH=$MPJ_HOME/bin:$JAVA_HOME/bin:$PATH```
	- Add library paths:```$ export LD_LIRARY_PATH=$MPJ_HOME/lib:$LD_LIBRARY_PATH ```
- Download the PanCGP package:
    - Extract it to your home directory, the path may be: _/home/< username>/PanCGP_
	- Now, we need to add path to environment variables
	- For PanCGP path, open .bashrc: ```$ vi .bashrc```
	- For adding path: ```$ export PANCGP=/home/<username>/PanCGP/depends/```
	- Add blast path: ```$ export PANCGPBLAST=/home/<username>/PanCGP/depends/blast/bin```
	- Add to path variable: ```$ export PATH=$PANCGP:$PANCGPBLAST:$PATH```
	- Set permissions of files in path to allow execution (with extension or without extension):
	    - ```$ sudo chmod +x *```
      - ```$ sudo chmod +x *.*```

    Now, PanCGP has been installed and ready to use for analysis purposes.

## User guide
PanCGP can be run in two modes: _Single-node_ and _multi-node_ cluster modes.
##### Single node
Assuming that user has successfully carried out the installation. To execute the PanCGP pipeline:
- As explained in installation, after extracting the PanCGP pipeline the path may be: _/home/< username>/PanCGP_
- Execute: ```$ cd /home/<username>/PanCGP```
- Execute: ```$ mpjrun.sh -np <no. of threads to run> -jar PanCGP.jar <Path to input folder containing FASTA files>```

##### Multi-node
This section list down the steps to run the PanCGP pipeline to run in a multi-node or cluster architecture.
- Step 1 and 2 from above section are same.
- To run PanCGP on multiple nodes, a _machine\_file_ is required, a text file listing all the IP or aliases of machines or nodes (which are fully qualified name of nodes in reality) on which pipeline is intended to execute.  <br />
    Suppose that there are three machines: machine1, machine2, machine3. In _machine\_file_, these will be listed as:

    machine1 <br />
    machine2 <br />
    machine3 <br />

- Execute: ```$ mpjboot <name of machine\_file>```  <br />
    This will start daemon processes on machines or nodes listed in the file. These daemon processes will host the parallel processes launched by executing PanCGP command.
- Execute: ```$ mpjrun.sh -np <total no. of threads to run> -dev niodev -jar PanCGP.jar <path to input folder containing FASTA files>```
- Once the pipeline has completed execution, execute: ```$ mpjhalt <name of machine\_file>```  <br />
    This command will stop the daemon processes running on machines or nodes listed in the file. There daemons are not necessary to stop after every execution. Once an execution is complete, they are ready to host another execution of pipeline.

### Version
1.0.0

### Todos

 - Write Tests
 - Rethink Github Save
 - Add Code Comments
 - Add Night Mode

License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [@thomasfuchs]: <http://twitter.com/thomasfuchs>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [marked]: <https://github.com/chjj/marked>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [keymaster.js]: <https://github.com/madrobby/keymaster>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]:  <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>


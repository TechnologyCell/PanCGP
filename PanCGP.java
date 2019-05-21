/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package PanCGP;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;
import mpi.*;
import org.apache.commons.io.FilenameUtils;
import java.security.MessageDigest;
import org.apache.commons.codec.binary.Hex;
import org.apache.commons.codec.binary.StringUtils;
/**
 *
 * @author ahsan
 */
public class PanCGP {

    private int Rank;
    private static final String BASEPATH = (PanCGP.class.getProtectionDomain().getCodeSource().getLocation().getPath()).substring(0, (PanCGP.class.getProtectionDomain().getCodeSource().getLocation().getPath()).lastIndexOf("/") + 1);
    private static final String DEPENDS = "";//BASEPATH + "depends";
    private static final String BLASTALL = DEPENDS + "blastall";
    private static final String TIGRCUT = DEPENDS + "tigrcut";
    private static final String FASTA_CONVERSION = DEPENDS + "saco_convert";
    private static final String GROUPING = DEPENDS + "Grouping.pl";
    private static final String NEWGENES = DEPENDS + "Newgenes.pl";
    private static final String NEWGENEFAMILIES = DEPENDS + "Newgenefamilies.pl";
    private static final String COREGENES = DEPENDS + "Coregenes.pl";
    private static final String GENEEXTRACTION = DEPENDS + "Subsets";

    private final String formatdb = DEPENDS + "formatdb";
    private final String cache_source = "blastp-core+matrix";

    private int totalThreads;
    
    private static long pid = 0;
    private String scratch;
    private String coregeneExtractionPath;
    private String pangenomeExtractionPath;
    private String dataPath;
    
    private int nprot = 0, skipped = 0, noOfFiles = 0, loopIterations = 0;
    private String checksumResult;
    
    private static Collection<HashMap<String, String>> config = new ArrayList<HashMap<String, String>>();
    private ArrayList<String> files = new ArrayList<String>();
    
    public PanCGP(ArrayList<String> filenameList, String dataPath, int mpiSize, int rank){
        this.files = filenameList;
        this.noOfFiles = this.files.size();
        if(dataPath.endsWith("/"))
            this.dataPath = dataPath;
        else
            this.dataPath = dataPath.concat("/");
        this.totalThreads = mpiSize;
        this.Rank = rank;
        
        this.loopIterations = (this.noOfFiles + this.totalThreads - 1)/this.totalThreads;
        
        String parentDir = this.dataPath.substring((this.dataPath.substring(0, this.dataPath.lastIndexOf("/")-1)).lastIndexOf("/")+1, this.dataPath.length()-1);
        scratch = this.dataPath + "Results/Pangenome."+ parentDir + "." + mpiSize + "/";
        coregeneExtractionPath = scratch + "coregenes/";
        pangenomeExtractionPath = scratch + "pangenomes/";
        
        File corePath = new File(coregeneExtractionPath);
        corePath.mkdirs();
        
        File panPath = new File(pangenomeExtractionPath);
        panPath.mkdirs();
        System.out.println("scratch = " + /*DEPENDS*/parentDir);
        //System.out.println("configsize = " + totalThreads);
    }
    public void parseConfig(){
        System.out.println("Parsing Fasta Files...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        for(int i=0; i<this.loopIterations; i++){
            HashMap<String,String> map = new HashMap<String, String>();
            int index = i * this.totalThreads + this.Rank;
            if (index < this.noOfFiles) {
                String name = files.get(index);
                map.put("description", name);
                map.put("source", name);
                map.put("id", Integer.toString(index));
                //System.out.println(index);
               // System.out.println(name);
                config.add(map);
            }
            //System.out.println(index);
        }
    }
    public void prepareFasta(){
        System.out.println("Hashing Data...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        String target;
        for(int i=0; i<this.loopIterations; i++){
            //Iterator<HashMap<String, String>> iterator = config.iterator();
            int index = i * this.totalThreads + this.Rank;
            if (index < this.noOfFiles && i < config.size()) {
                String tempWriter = scratch + "checksum("+Integer.toString(index)+").totalgenes.list";
                try{    
                    File filetoWrite = new File(tempWriter);
                    if (!filetoWrite.exists()){
                        filetoWrite.createNewFile();
                    }
                    FileWriter writer;
                    writer = new FileWriter(filetoWrite);
                    BufferedWriter bw = new BufferedWriter(writer);

                    fasta2Fasta(index);
                    HashMap<String, String> map = (HashMap<String, String>) config.toArray()[i];

                    target = scratch + index + ".fsa";
                    map.put("target", target);
                    map.put("total genes", Integer.toString(nprot));
                    map.put("skipped", Integer.toString(skipped));
                    map.put("checksum", checksumResult);
//                        System.out.println("Target => " + target);
//                        System.out.println("Total Genes => " + nprot);
//                        System.out.println("Skipped => " + skipped);
//                        System.out.println("Hash Result => " + checksumResult);

                    bw.write(checksumResult + "." + Integer.toString(nprot) + "\n");
                    if(nprot == 0){
                        //////  Exit the scenerio  ////////////
//                            System.out.println("|||||||||||||| NO GENES FOUND ||||||||||||||||||||||||||");
                    }
                    File file = new File(scratch + "group_" + Integer.toString(index) + ".dat");
                    if (!file.exists()) file.createNewFile();
                    bw.close();
                } catch (IOException ex) {
                    Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    } 
    public void fasta2Fasta(int id){
        String tempReader = scratch + id + ".fsa.temp.raw";
        String tempWriter = scratch + id + ".fsa";
        String filename = dataPath + files.get(id);
        
        //Initializing for every Strain
        nprot = skipped = 0;
        
        try{
            String[] cmd = {"/bin/sh",
                            "-c",
                            FASTA_CONVERSION + " -I fasta -O raw " + filename + " | sort > "+tempReader
                            };
            System.out.println(cmd[2]);
            Process p = Runtime.getRuntime().exec(cmd);
            int r = p.waitFor();
//            System.out.println("result = " + r);
            //File file = new File(dataPath+files.get(id));
            ////////////////////////////////////////// MD5 Hash  ////////////////////////////////////////
            File file = new File(tempReader);
            ArrayList<String> filelist = new ArrayList<String>();
            filelist.add(tempReader);
            md5(filelist);
//            System.out.println(checksumResult);
            ////////////////////////////////////////////////////////////////////////////////////////////
            
            File filetoWrite = new File(tempWriter);
            if (!file.exists()){
                file.createNewFile();
            }
            FileWriter writer = new FileWriter(filetoWrite);
            BufferedWriter bw = new BufferedWriter(writer);
            File filetoRead = new File(tempReader);
            BufferedReader fileReader = new BufferedReader(new FileReader(filetoRead));
            String tmpstr = null;
            while((tmpstr = fileReader.readLine()) != null){
                if(!tmpstr.matches("([A-Za-z]+)")){
                    skipped++;
                    continue;
                }
                nprot++;
                bw.write(">"+ checksumResult +"."+ nprot +"\n");
                for(int i=0; i<tmpstr.length(); i+=60) {
                    if(tmpstr.length() < i + 60)
                        bw.write(tmpstr.substring(i, i+(tmpstr.length()%60))+"\n");
                    else
                        bw.write(tmpstr.substring(i, i+60) +"\n");
                }
                //bw.write(tmpstr+"\n");
            }
            bw.close();
            
        }catch(Exception e)
        {
            System.out.println(e.toString());
        }
        
       // System.out.println("/tmp/Pangenome.13215.0/0.fsa.tempReader.raw");
    }
    public String md5(ArrayList<String> files){
        try{
            int i = 0;
            String tmpstr = null, line = null;
            MessageDigest digest = MessageDigest.getInstance("MD5");
            for(int j=0; j<files.size(); j++){
                BufferedReader reader = new BufferedReader(new FileReader(files.get(j)));
                while((tmpstr = reader.readLine()) != null){
                    if (i == 0) {
                        line = tmpstr;
                        i++;
                    }
                    line.concat(tmpstr);
                    digest.update(tmpstr.getBytes());
                   //System.out.println(line);
                }
            }
            byte[] resultByte = digest.digest();
            checksumResult = new String(Hex.encodeHex(resultByte));
        }catch(Exception e){
            System.out.println(e.toString());
        }
        return checksumResult;
    }
    
    public void makeBlastDB() {
        System.out.println("Creating BLAST Databases...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        int i = 0, index = 0;
       // System.out.println(config.size() + " size with rank " + this.Rank);
        try{
            for(i=0; i<this.loopIterations; i++){
//                Rank = i;
                index = i * this.totalThreads + this.Rank;
                if (index < this.noOfFiles && i < config.size()) {
                    HashMap<String, String> map = (HashMap<String, String>) config.toArray()[i];
                    String target = map.get("target");
                    //System.out.println(target);
                    String[] cmd = {"/bin/sh",
                            "-c",
                            formatdb + " -i " + target + " -p T -t " + Integer.toString(index)
                            };

                    System.out.println(cmd[2]);
                    Process p = Runtime.getRuntime().exec(cmd);
                    int a = p.waitFor();
                    //System.out.println(a + " " + i + " with rank " + index + " " + map.get("source"));
                }
            }
        }catch(Exception e){
            System.out.println(e.toString() + " at iteration => " + Integer.toString(index));
        }
    } 
    public void makeBlastReports() throws IOException{
        System.out.println("Running BLAST...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        int index = 0;
        for(int i=0; i<this.loopIterations; i++){
            for(int j=0; j<this.noOfFiles; j++){
                index = i * this.totalThreads + this.Rank;
                if (index < this.noOfFiles) {
                    String a = scratch + Integer.toString(index) + ".fsa";
                    String b = scratch + Integer.toString(j) + ".fsa";
                    String blast_scratch = scratch + Integer.toString(index) + "-" + Integer.toString(j) + ".blastout.gz";
                    ArrayList<String> filelist = new ArrayList<String>();
                    filelist.add(a);
                    filelist.add(b);
                    String jobId = md5(filelist);
                    System.out.println("Job id = " + jobId);
                     String[] cmd = {"/bin/sh",
                            "-c",
                            BLASTALL + " -F 0 -i " + a + " -p blastp -e 1e-5 -m 7 -d " + b + " | " + TIGRCUT + " | gzip > " + blast_scratch
                            };
                    System.out.println(cmd[2]);
                    Process p = Runtime.getRuntime().exec(cmd);
                    try {
                        int ab = p.waitFor();
                        //System.out.println(ab);
                    } catch (Exception ex) {
                        Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        }
    }
    
    public void preGroup() {
        System.out.println("Preparing Hashing For Clustering...\n");
        String tempWriter = scratch + "checksum.totalgenes.list";
        try{    
            File filetoWrite = new File(tempWriter);
            if (!filetoWrite.exists()){
                filetoWrite.createNewFile();
            }
            FileWriter writer;
            writer = new FileWriter(filetoWrite);
            BufferedWriter bw = new BufferedWriter(writer);
            
            String tempReader;
            for(int i=0; i<this.noOfFiles; i++){
                tempReader = scratch + "checksum("+Integer.toString(i) +").totalgenes.list";
                File filetoRead = new File(tempReader);
                if (!filetoRead.exists()){
                    System.out.println("File not found");
                }
                FileReader reader;
                reader = new FileReader(filetoRead);
                BufferedReader br = new BufferedReader(reader);
                String tmpStr = br.readLine();
                
                bw.write(tmpStr+"\n");
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void clustering() {
        //String scratch = "/tmp/Pangenome.2257.0";
        System.out.println("Clustering...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        int index = 0;
        for(int i=0; i<this.loopIterations; i++) {
            index = i * this.totalThreads + this.Rank;
            if (index < this.noOfFiles) {
                try{
                    String[] cmd = {"/bin/sh",
                        "-c",
                        GROUPING + " -dirpath " + scratch + " -rank " + Integer.toString(index) + " " + scratch + "checksum.totalgenes.list"
                        };
                    System.out.println(cmd[2]);
                    Process p = Runtime.getRuntime().exec(cmd);
                    //int a = p.waitFor();
                   // System.out.println(a);
                }
                catch (Exception ex) {
                    Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }
    
    public void panGenome() {
        System.out.println("Calculating PanGenome...\n");
        String tempReader, filepangenomes =  scratch + "pangenomes.txt";
        
        try{ 
            File file = new File(filepangenomes);
            if (!file.exists())     file.createNewFile();
            
            FileWriter writer;
            writer = new FileWriter(file);
            BufferedWriter bw = new BufferedWriter(writer);

            for(int i=0; i<this.noOfFiles; i++) {
                int pangenome = 0;
                tempReader = scratch + "group_" + Integer.toString(i) + ".dat";
                File filetoRead = new File(tempReader);
                if (!filetoRead.exists()){
                    System.out.println("File not found");
                }
                FileReader reader;
                reader = new FileReader(filetoRead);
                BufferedReader br = new BufferedReader(reader);

                while((br.readLine()) != null)  pangenome++;
                
                bw.write(Integer.toString(pangenome) + "\n");
                br.close();
                //System.out.println(pangenome);
            } 
            bw.close();
        } catch (Exception ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void newGenes() {
        System.out.println("Calculating New Genes...\n");
        Iterator<HashMap<String, String>> iterator = config.iterator();
        String tempReader, tmpStr, filenewgenes;
        filenewgenes = scratch + "newgenes.txt";
        
        try{    
            File file = new File(filenewgenes);
            if (!file.exists())     file.createNewFile();
            
            String[] cmd = {"/bin/sh",
                "-c",
                NEWGENES + " -dirpath " + scratch + " " + scratch + "checksum.totalgenes.list"
                };
            System.out.println(cmd[2]);
            Process p = Runtime.getRuntime().exec(cmd);
            int a = p.waitFor();
            //System.out.println(a);
            
            tempReader = filenewgenes;
            File filetoRead = new File(tempReader);
            FileReader reader;
            reader = new FileReader(filetoRead);
            BufferedReader br = new BufferedReader(reader);
            br.close();
        } catch (Exception ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void coreGenes() {
        System.out.println("Calculating Core Genes...\n");
        Iterator<HashMap<String, String>> iterator = config.iterator();
        String tempReader, tmpStr, filecoregenes;
        filecoregenes = scratch + "coregenes.txt";
        
        try{
            File file = new File(filecoregenes);
            if (!file.exists())     file.createNewFile();
            
            String[] cmd = {"/bin/sh",
                "-c",
                /*"perl " +*/ COREGENES + " -dirpath " + scratch + " " + scratch + "checksum.totalgenes.list"
                };
            System.out.println(cmd[2]);
            Process p = Runtime.getRuntime().exec(cmd); 
            int a = p.waitFor();
            //System.out.println(a);
            
            tempReader = filecoregenes;
            File filetoRead = new File(tempReader);
            FileReader reader;
            reader = new FileReader(filetoRead);
            BufferedReader br = new BufferedReader(reader);
            br.close();
        } catch (Exception ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }        
    }
    public void newGeneFamilies() {
        System.out.println("Calculating New Gene Families...\n");
        String tempReader, tmpStr, filenewgenefamilies;
        filenewgenefamilies = scratch + "newgenefamilies.txt";

        try{    
            File file = new File(filenewgenefamilies);
            if (!file.exists())     file.createNewFile();
            
            String[] cmd = {"/bin/sh",
                "-c",
                NEWGENEFAMILIES + " -dirpath " + scratch + " " + scratch + "checksum.totalgenes.list"
                };
            System.out.println(cmd[2]);
            Process p = Runtime.getRuntime().exec(cmd);
            int a = p.waitFor();
            //System.out.println(a);
            
            tempReader = filenewgenefamilies;
            File filetoRead = new File(tempReader);
            FileReader reader;
            reader = new FileReader(filetoRead);
            BufferedReader br = new BufferedReader(reader);
            br.close();
        } catch (Exception ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }                
    }
    public void makeTable() {
        System.out.println("Compiling Results...\n");
        try{
            String readpan = scratch + "pangenomes.txt",
                   readcore = scratch + "coregenes.txt",
                   readnewgenes = scratch + "newgenes.txt",
                   readfamilies = scratch + "newgenefamilies.txt",
                   readtotal = scratch + "checksum.totalgenes.list",
                    
                   writeTable = scratch + "table.txt";
            
            File filepan = new File(readpan), 
                 filecore = new File(readcore), 
                 filenewgenes = new File(readnewgenes),
                 filefamilies = new File(readfamilies),
                 filetotal = new File(readtotal),
                    
                 fileTable = new File(writeTable);
            
            FileReader readerpan = new FileReader(filepan),
                       readercore = new FileReader(filecore),
                       readernewgenes = new FileReader(filenewgenes),
                       readerfamilies = new FileReader(filefamilies),
                       readertotal = new FileReader(filetotal);
            
            BufferedReader brpan = new BufferedReader(readerpan),
                           brcore = new BufferedReader(readercore),
                           brnewgenes = new BufferedReader(readernewgenes),
                           brfamilies = new BufferedReader(readerfamilies),
                           brtotal = new BufferedReader(readertotal);
            
            FileWriter writertable= new FileWriter(writeTable);
            BufferedWriter bwtable = new BufferedWriter(writertable);
            
            if(!fileTable.exists()) fileTable.createNewFile();
            
            bwtable.write("Prot.No." + "\t" + "Names" + "\t" + "Total Proteins" + "\t" + "Unique gene/proteins" + "\t" + "New gene/protein Clusters" + "\t" + "Pangenome" + "\t" + "Core genome/proteome" + "\n");
            
            for(int i=0; i<this.noOfFiles; i++) {
                String str = files.get(i);
            }
            
            for(int i=0; i<this.noOfFiles; i++) {
                String filename = FilenameUtils.getBaseName(FilenameUtils.getBaseName(files.get(i))).replaceAll("_prodigal", "");
                String[] filenamenew = filename.split("_");
                String finalFileName = "";
                int j=0;
                while(j< filenamenew.length){
                    String str = filenamenew[j];
                    
                    if(str.matches("^[a-zA-Z]+$") && !str.matches("^[0-9]+")) {
                        finalFileName += str.substring(0, 1);
                        System.out.println(str);
                    }else{
                        finalFileName += "_" + str;
                    }
                    //System.out.println("filename = > " + finalFileName);
                    j++;
                }
                bwtable.write(Integer.toString(i) + "\t" + finalFileName + "\t" + brtotal.readLine().split("\\.")[1] + "\t" + brnewgenes.readLine() + "\t" + brfamilies.readLine() + "\t" + brpan.readLine() + "\t" + brcore.readLine() + "\n");
            }
            bwtable.close();
            
            brpan.close();
            brcore.close();
            brnewgenes.close();
            brfamilies.close();
            brtotal.close();
        } catch (Exception ex) {
            Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("Successful\nResults are placed in " +  scratch + "table.txt");
    }
    
    public void coreGenesExtraction() {
        System.out.println("Core Gene Extraction...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        String filecoregene;
        int index = 0;
        for(int i=0; i<this.loopIterations; i++) {
            index = i * this.totalThreads + this.Rank;
            if (index < this.noOfFiles) {
                try{   
                    filecoregene = coregeneExtractionPath + this.files.get(index);
                    String[] cmd = {"/bin/sh",
                        "-c",
                        /*"perl " +*/ GENEEXTRACTION + " -i 1:" + Integer.toString(index+1) + " " + scratch + " > " + filecoregene
                        };
                    System.out.println(cmd[2]);
                    Process p = Runtime.getRuntime().exec(cmd);
                    int a = p.waitFor();
                    //System.out.println(a);
                } catch (Exception ex) {
                    Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }
    
    public void panGenomeExtraction() {
        System.out.println("PanGenome Extraction...Invoking Thread " + Integer.toString(this.Rank) + "\n");
        String filepangenome;
        int index = 0;
        for(int i=0; i<this.loopIterations; i++) {
            index = i * this.totalThreads + this.Rank;
            if (index < this.noOfFiles) {
                try{   
                    filepangenome = pangenomeExtractionPath + this.files.get(index);
                    String[] cmd = {"/bin/sh",
                        "-c",
                        GENEEXTRACTION + " -u 1:" + Integer.toString(index+1) + " " + scratch + " > " + filepangenome
                        };
                    System.out.println(cmd[2]);
                    Process p = Runtime.getRuntime().exec(cmd);
                    int a = p.waitFor();
                    //System.out.println(a);
                } catch (Exception ex) {
                    Logger.getLogger(PanCGP.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
    }
    
    public void getCollection(String path) {
        dataPath = path;
        File actual = new File(path);
        for( File f : actual.listFiles()){
            String ext = FilenameUtils.getExtension(f.getName());
            if(ext.contains("fsa")){
                files.add(f.getName());
               // System.out.println( f.getName()+" with "+ext );
            }
        }
    }
    
    public static void main(String[] args) throws IOException, InterruptedException {

        long starttime = 0;
        MPI.Init(args);
        int me = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();
  //      System.out.println("hi from < " + me +" >");
        if (me == 0) {
            System.out.println("\nRunning PanCGP on " + size + " threads... \n");
            starttime = System.currentTimeMillis();
        }
        
        
        ArrayList<String> files = new ArrayList<String>();
        File actual = new File(args[3]);
        for( File f : actual.listFiles()){
            String ext = FilenameUtils.getExtension(f.getName());
            if(ext.contains("fsa")){
                files.add(f.getName());
               // System.out.println( f.getName()+" with "+ext );
            }
        }
        PanCGP test = new PanCGP(files, args[3], size, me);

        MPI.COMM_WORLD.Barrier();
        test.parseConfig();
        
        MPI.COMM_WORLD.Barrier();
 
        test.prepareFasta();
        MPI.COMM_WORLD.Barrier();
        
        test.makeBlastDB();
        MPI.COMM_WORLD.Barrier();
        
        test.makeBlastReports();
        if(me == 0)
            test.preGroup();
        MPI.COMM_WORLD.Barrier();
        
        test.clustering();   
        MPI.COMM_WORLD.Barrier();

        if(me == 0) {
            test.panGenome();
            test.coreGenes();
            test.newGeneFamilies();
            test.newGenes();
            test.makeTable();
        }
        MPI.COMM_WORLD.Barrier();
        
        test.coreGenesExtraction();
        test.panGenomeExtraction();
        MPI.Finalize(); 
        if(me == 0) {
            System.out.println(Integer.toString((int)((System.currentTimeMillis() - starttime)/1000)) + " secs Elapsed");
        }
    }
}

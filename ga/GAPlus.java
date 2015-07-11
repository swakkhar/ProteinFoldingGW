/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ga;

/**
 *
 * @author Administrator
 */
import java.util.ArrayList;
import java.util.Collections;
import java.io.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class GAPlus {
    private AminoAcid aa3d;
    private Conformation c3d;    
    int seqLen = 0;
    int latticeModel;
    boolean needRestart=false;
    Global g = new Global();
    private ArrayList<Conformation> massPop=new ArrayList<Conformation>();
    private String sOutput="";
    Operators _operator;
    int _pctAllowautoTryCount=100;
    public GAPlus(String targetSeq) throws Exception {
        
        _operator=new Operators();
        aa3d = new AminoAcid();
        c3d = new Conformation();        
        Global.seqLength = targetSeq.length();
        Global.aaseq=targetSeq;        
        
        seqLen = Global.seqLength;
        Global.startTimeInSec=System.currentTimeMillis()/1000;
        
        //Global.arraySeq=Global.inputSequence_01.toCharArray();
        //xy-plane @z=0 and denoted the vector as a,b,c, and d
        Global.vector3Dfcc[0] = new Node3d(+1, +1, +0);
        Global.vector3Dfcc[1] = new Node3d(+1, -1, +0);
        Global.vector3Dfcc[2] = new Node3d(-1, -1, +0);
        Global.vector3Dfcc[3] = new Node3d(-1, +1, +0);

        //yz-plane @x=0 and denoted the vectors as e,f,g, and h
        Global.vector3Dfcc[4] = new Node3d(+0, +1, +1);
        Global.vector3Dfcc[5] = new Node3d(+0, -1, +1);
        Global.vector3Dfcc[6] = new Node3d(+0, -1, -1);
        Global.vector3Dfcc[7] = new Node3d(+0, +1, -1);
        //zx-plane @y=0 and denoted the vectors as i,j,k, and l
        Global.vector3Dfcc[8] = new Node3d(+1, +0, +1);
        Global.vector3Dfcc[9] = new Node3d(+1, +0, -1);
        Global.vector3Dfcc[10] = new Node3d(-1, +0, -1);
        Global.vector3Dfcc[11] = new Node3d(-1, +0, +1);

        //12 vector array
        //ArrayList<ArrayList> myArrays=new ArrayList<ArrayList>();
        Global.vecArrays.clear();
        Global.vecArrays = FCCDescriptor.get12x12BasisVectors();
        Global.vecA = Global.vecArrays.get(0);
        Global.vecB = Global.vecArrays.get(1);
        Global.vecC = Global.vecArrays.get(2);
        Global.vecD = Global.vecArrays.get(3);
        Global.vecE = Global.vecArrays.get(4);
        Global.vecF = Global.vecArrays.get(5);
        Global.vecG = Global.vecArrays.get(6);
        Global.vecH = Global.vecArrays.get(7);
        Global.vecI = Global.vecArrays.get(8);
        Global.vecJ = Global.vecArrays.get(9);
        Global.vecK = Global.vecArrays.get(10);
        Global.vecL = Global.vecArrays.get(11);
    }

public void startSearch(int pctAllowSimilar) throws Exception{
    _pctAllowautoTryCount=pctAllowSimilar;
    int amdschk=0;
    sOutput="--rndSeed: "+Global.rndSeed+" popSize: "+Global.populationSize+", seqLen: "+Global.seqLength+", pctSimilar: "+pctAllowSimilar+" *********************";
    new PrintStream(Global.fout).println (sOutput);    

    if (Global.traceOn) System.out.println(sOutput);
    sOutput="--***********************************************************************************************";
    new PrintStream(Global.fout).println (sOutput);    

    if (Global.traceOn) System.out.println(sOutput);
    
    AminoAcid[] aa=new AminoAcid[Global.seqLength];
    ArrayList<Conformation> newPop=new ArrayList<Conformation>();
    Conformation c;
    int lastImpGen=0;
    long nonImpRestartInSec=0;
    double primaryf;
    int secondaryf;
    Global.generationCount=0;
    Global.restartCount=0;
//    Global.startTimeInSec=0;
    long currTimeInSec=0;
    long stopTimeInSec=0;
    long elapseTimeSec=0;
    long ioTimeMilliSec=0;
    long ioStartMilliSec=0;
    long ioEndMilliSec=0;
    long improveTimeInSec=0;
    int completeRestartCount=0;
    int similar=0;
    int similarLimit=0;
    int shuffleincr=0;
    int avgtracking=0;

    double avf=0,tpf=0;
    double gbf=0;
    int cbf=0;
    int cbfsim=0;
    int counter=0;
    boolean isCompleteRestart=false;
    String enSeq="";
    needRestart=false;
//    _pctAllowSimilar=pctAllowSimilar;
    Global.generationCount=0;
    Global.bestFitConf=Global.arrayCurrentPopulation[0];
    Global.bestFitness=Global.bestFitConf.energy;
    primaryf=Global.bestFitConf.energy;
    gbf=primaryf;

    secondaryf=Global.bestFitConf.fitness[1];
    cbf=Global.bestFitConf.fitness[0];
    cbfsim=Global.bestFitConf.fitness[1];
    enSeq=Global.bestFitConf.encodedString;


//    Global.startTimeInSec=System.currentTimeMillis()/1000;
    try{
        sOutput="('"+Global.pid+"',"+lastImpGen+","+gbf+","+elapseTimeSec+","+Global.restartCount+","+Global.solExplored+",'"+enSeq+"',"+Global.discardedTwins+");";
        new PrintStream(Global.fout).println (sOutput);
            
        shuffleincr=0;
        similarLimit=Global.seqstagparam;
        //int greedtracking=0;
        
        //op logging
        int opusedlast=-1;
        double improvment=0;        
        
        while (counter<=Global.maxTry){            
            tpf=0;
            for (int p=0;p<Global.populationSize;p++){
                tpf=tpf+Global.arrayCurrentPopulation[p].energy;
            }
            avf=tpf/Global.populationSize;
            counter++;
            currTimeInSec=System.currentTimeMillis()/1000;
            int diff=0;
            int ra=Common.getRandomNext(0, 100);
            if(ra<25){ //crossover
                if (ra<5){//10                    
                    crossoverGenerationMulti();
                    opusedlast=0;
                }
                else{//20                    
                    crossoverGeneration();
                    opusedlast=1;
                }                
            }
            else{ //mutation                                
                if (ra<40){//30
                    mutatedGeneration();
                    opusedlast=2;
                }
                else if(ra<55){//40
                    pullSinglePointGeneration();
                    opusedlast=3;
                }
                else if(ra<75){//50
                    cornerplipGeneration();
                    opusedlast=4;
                }
                else if(ra<100){//70
                    pullGeneration();
                    opusedlast=5;
                }                
                /* else if(ra<90){//90
                    tiltGeneration();
                    opusedlast=6;
                 }                                
                else if(ra<100){
                    for (int p=1;p<Global.populationSize;p++){
                        if(Global.arrayCurrentPopulation[p].crankshaft)
                            continue;
                        Global.arrayCurrentPopulation[p].crankshaft=true;

                        Conformation cx=new Conformation(Global.arrayCurrentPopulation[p]);
                        cx=_operator.crankShaftRotation(cx, 4);
                        if(cx.compareTo(Global.arrayCurrentPopulation[p])<0)
                            Global.arrayCurrentPopulation[p]=new Conformation(cx);
                    }
                    Common.sortByFitness(Global.arrayCurrentPopulation);

                    opusedlast=7;
                }*/
            }            

            if (cbfsim==Global.arrayCurrentPopulation[0].fitness[1]){
                similar++;
            }else{
                cbfsim=Global.arrayCurrentPopulation[0].fitness[1];
                similar=0;
            }
            cbf=Global.arrayCurrentPopulation[0].fitness[0];
            
            
            if (cbf<Global.bestFitConf.fitness[0]){
                //simulated fixed but originally improving let improve
                 similar=0;
                 improvment=Global.arrayCurrentPopulation[0].energy-Global.bestFitness;

                 shuffleincr=0;

                 Global.bestFitConf=Global.arrayCurrentPopulation[0];
                 primaryf=Global.bestFitConf.energy;
                 secondaryf=Global.bestFitConf.fitness[1];
                 
                 Global.bestFitness=primaryf;

                 lastImpGen=counter;
                 enSeq=Global.bestFitConf.encodedString;

                 currTimeInSec=System.currentTimeMillis()/1000;
                 elapseTimeSec=currTimeInSec-Global.startTimeInSec;
                 improveTimeInSec=elapseTimeSec;
                 nonImpRestartInSec=elapseTimeSec;
                 isCompleteRestart=false;

                 ioStartMilliSec=System.currentTimeMillis();
                 
                 gbf=primaryf;
                 
                 sOutput="('"+Global.pid+"',"+lastImpGen+","+gbf+","+elapseTimeSec+","+Global.restartCount+","+Global.solExplored+",'enSeq',"+Global.discardedTwins+");";
                 new PrintStream(Global.fout).println (sOutput);
                 ioEndMilliSec=System.currentTimeMillis();
                 ioTimeMilliSec+=ioEndMilliSec-ioStartMilliSec;
                // if (Global.traceOn) System.out.println(sOutput);
             }else{
                improvment=0;
             }

             currTimeInSec=System.currentTimeMillis()/1000;
             elapseTimeSec=currTimeInSec-Global.startTimeInSec;

             if (Global.traceOn||Global.logging){
                 double tmpavf=avf*1000;
                 tmpavf=Math.round(tmpavf);
                 tmpavf/=1000;
                  sOutput="('"+Global.pid+"',"+counter+","+primaryf+","+cbf+","+tmpavf+","+secondaryf+","+ Global.restartCount+","+elapseTimeSec+","+ioTimeMilliSec+","+Global.solExplored+","+Global.discardedTwins+");";
                 if(Global.traceOn) System.out.println(sOutput+" "+1);
                 if(Global.logging) new PrintStream(Global.foutlog).println (sOutput);
                 
                  //op logging
                  sOutput="("+counter+","+opusedlast+","+improvment+")";
                 if(Global.traceOn) System.out.println(sOutput+" "+1);
                 if(Global.logging) new PrintStream(Global.foutoplog).println (sOutput);                 
            }            
             
            //===========stagnation recovery=================================================
            if(similar>similarLimit){
                /*if(greedtracking>Global.distswitchparam){
                    switchEnergyDist();
                    greedtracking=0;
                    shuffleincr=0;                    
                }
                else{*/
                    if(shuffleincr<3)
                        shuffleincr++;
                    else
                        shuffleincr=1;
                    
                    sOutput="Stagnation recovery: "+shuffleincr;
                    if (Global.traceOn){
                        System.out.println(sOutput);
                    }
                    if(Global.logging){
                        new PrintStream(Global.foutlog).println (sOutput);
                    }

                    improveTimeInSec=elapseTimeSec;
                    lastImpGen=counter;
                    Global.restartCount+=1;
                    int shuffling=5*shuffleincr;
                    shuffling+=Global.ranDOM.nextInt(shuffling);
                    go4randomWalk(shuffling);
                //}
            }

            isCompleteRestart=false;
            stopTimeInSec=System.currentTimeMillis()/1000;
            diff=(int)(stopTimeInSec-Global.startTimeInSec);

            //param for dist switch
            //greedtracking++;
        }
        sOutput="('"+Global.pid+"',"+counter+","+gbf+","+elapseTimeSec+","+Global.restartCount+","+Global.solExplored+",'"+enSeq+","+Global.discardedTwins+"');";

        new PrintStream(Global.fout).println (sOutput);

        new PrintStream(Global.foutlog).println (sOutput);
        new PrintStream(Global.foutlog).println (sOutput+" "+secondaryf);
        
        if (Global.traceOn){
            System.out.println(sOutput);
            System.out.println(sOutput+" "+secondaryf);
        }

        for(int i=0; i<Global.populationSize; i++){
            sOutput=Global.arrayCurrentPopulation[i].encodedString;
            new PrintStream(Global.fout).println (sOutput);
        }
    }
    catch (Exception e) {
        //System.out.println(e.getMessage());
        System.out.println("Error @startSearch: "+e.toString());
    }
}
private void switchEnergyDist(){
    Global.energydistselected=1-Global.energydistselected;
    
    for(int i=0; i<Global.populationSize; i++){
        Common.calculateFitness(Global.arrayCurrentPopulation[i].aminoAcid,Global.arrayCurrentPopulation[i].fitness);
        Global.arrayCurrentPopulation[i]=new Conformation(Global.arrayCurrentPopulation[i].aminoAcid, Global.arrayCurrentPopulation[i].fitness, Global.arrayCurrentPopulation[i].encodedString);
    }
    sOutput="GreedMode"+Global.energydistselected;
    if(Global.traceOn) System.out.println(sOutput+" "+1);
    if(Global.logging) new PrintStream(Global.foutlog).println (sOutput);
}

private void go4randomWalk(int _shuffling) throws Exception{
    for (int k=0;k<_shuffling;k++){
        int j=0;

        for(Conformation el:Global.arrayCurrentPopulation){
            for (int pos=Common.getRandomNext(0,Global.seqLength-3)+1;pos<Global.seqLength;pos++)
            {
                Conformation cx2=_operator.randomWalk(el, pos);
                if (cx2==null) continue;
                Global.arrayCurrentPopulation[j]=new Conformation(cx2);
                break;
            }
            j++;
        }
    }
}
private void tiltGeneration() throws Exception{
    int pos;
    Conformation c=new Conformation();
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    for(int p=0;p<Global.populationSize;p++){
        neighbor.clear();
         if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }
        trycount++;
        neighbor.add(Global.arrayCurrentPopulation[p]);

        if(Global.arrayCurrentPopulation[p].tilt){
            pop.add(Global.arrayCurrentPopulation[p]);
            continue;
        }
        Global.arrayCurrentPopulation[p].tilt=true;

        for (int i=1;i<Global.seqLength-3;i++){
            c=_operator.tiltMove(Global.arrayCurrentPopulation[p], i,i+1);
            if (c!=null) {
                neighbor.add(new Conformation(c));
            }
          }
        Collections.sort(neighbor);

        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
    }
    Collections.sort(pop);
    int i=0;
    for (Conformation el: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(el);
        if (i>=Global.populationSize) break;
    }
}
private void crossoverGenerationMulti() throws Exception{
    int idx1,idx2;
    boolean isExist=false;
    Conformation[] cout=new Conformation[2];
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    while(pop.size()<Global.populationSize){
        neighbor.clear();
        if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }
        trycount++;
        idx1=Common.getRandomNext(0, Global.populationSize-1);
        idx2=Common.getRandomNext(0, Global.populationSize-1);
        neighbor.add(Global.arrayCurrentPopulation[idx1]);
        neighbor.add(Global.arrayCurrentPopulation[idx2]);
        for (int i=1;i<Global.seqLength/10;i++){
            //int pos=Common.getRandomNext(2, Global.seqLength-2);
            _operator.crossOverMulti(idx1,idx2, cout);
            //if (cout==null) continue;
            if (cout[0]!=null) {
                neighbor.add(new Conformation(cout[0]));
            }
            if (cout[1]!=null) {
                neighbor.add(new Conformation(cout[1]));
            }
        }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
        if (!Common.isExist(neighbor.get(1),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(1));
    }
    Collections.sort(pop);
//    System.out.print("xover: "+pop.get(0).fitness);
    int i=0;
    for (Conformation c: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(c);
        if (i>=Global.populationSize) break;
    }
}
private void crossoverGeneration() throws Exception{
    int idx1,idx2;
    boolean isExist=false;
    Conformation[] cout=new Conformation[2];
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    while(pop.size()<Global.populationSize){
        neighbor.clear();
        if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }
        trycount++;
        idx1=Common.getRandomNext(0, Global.populationSize-1);
        idx2=Common.getRandomNext(0, Global.populationSize-1);
        neighbor.add(Global.arrayCurrentPopulation[idx1]);
        neighbor.add(Global.arrayCurrentPopulation[idx2]);
        for (int i=1;i<Global.seqLength/10;i++){
            int pos=Common.getRandomNext(2, Global.seqLength-2);
            _operator.crossOver(idx1,idx2, pos, cout);
            //if (cout==null) continue;
            if (cout[0]!=null) {
                neighbor.add(new Conformation(cout[0]));
            }
            if (cout[1]!=null) {
                neighbor.add(new Conformation(cout[1]));
            }
        }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
        if (!Common.isExist(neighbor.get(1),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(1));
    }
    Collections.sort(pop);
    int i=0;
    for (Conformation c: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(c);
        if (i>=Global.populationSize) break;
    }
}

private void cornerplipGeneration() throws Exception{
    int pos;
    boolean isExist=false;
    Conformation c=new Conformation();
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    for(int p=0;p<Global.populationSize;p++){
        neighbor.clear();
         if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }
        trycount++;
        neighbor.add(Global.arrayCurrentPopulation[p]);

        if(Global.arrayCurrentPopulation[p].cornerflip){
            pop.add(Global.arrayCurrentPopulation[p]);
            continue;
        }
        Global.arrayCurrentPopulation[p].cornerflip=true;

        for (int i=0;i<Global.seqLength;i++){
            c=_operator.cornerFlip(Global.arrayCurrentPopulation[p],i);
            if (c!=null) {
                neighbor.add(new Conformation(c));
            }
          }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
    }
    Collections.sort(pop);
    int i=0;
    for (Conformation el: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(el);
        if (i>=Global.populationSize) break;
    }
}

private void pullGeneration() throws Exception{
    int pos;
    Conformation c=new Conformation();
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    for(int p=0;p<Global.populationSize;p++){
        neighbor.clear();
         if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }

        trycount++;

        neighbor.add(Global.arrayCurrentPopulation[p]);

        if(Global.arrayCurrentPopulation[p].pull){
            pop.add(Global.arrayCurrentPopulation[p]);
            continue;
        }
        Global.arrayCurrentPopulation[p].pull=true;

        for (int i=0;i<Global.seqLength;i++){
            c=_operator.pullMove(Global.arrayCurrentPopulation[p], i);
            if (c!=null) {
                neighbor.add(new Conformation(c));
            }
          }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
    }

    Collections.sort(pop);
    int i=0;
    for (Conformation el: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(el);
        if (i>=Global.populationSize) break;
    }
}
private void pullSinglePointGeneration() throws Exception{
    int pos;
    Conformation c=new Conformation();
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    for(int p=0;p<Global.populationSize;p++){
        neighbor.clear();
         if(trycount>=Global.populationSize*10){
           // needRestart=true;
            break;
        }
        trycount++;
        neighbor.add(Global.arrayCurrentPopulation[p]);

        if(Global.arrayCurrentPopulation[p].pullsinglepoint){
            pop.add(Global.arrayCurrentPopulation[p]);
            continue;
        }
        Global.arrayCurrentPopulation[p].pullsinglepoint=true;

        for (int i=0;i<Global.seqLength;i++){
            c=_operator.pullMoveSinglePoint(Global.arrayCurrentPopulation[p],i);
            if (c!=null) {
                neighbor.add(new Conformation(c));
            }
          }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
    }
    Collections.sort(pop);
    int i=0;
    for (Conformation el: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(el);
        if (i>=Global.populationSize) break;
    }
}
private void mutatedGeneration() throws Exception{
    int pos;
    Conformation c=new Conformation();
    ArrayList<Conformation> pop=new ArrayList<Conformation>();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    int trycount=0;
    for(int p=0;p<Global.populationSize;p++){
        neighbor.clear();
         if(trycount>=Global.populationSize*10){
            //needRestart=true;
            break;
        }
        trycount++;
        neighbor.add(Global.arrayCurrentPopulation[p]);

        if(Global.arrayCurrentPopulation[p].mutate){
            pop.add(Global.arrayCurrentPopulation[p]);
            continue;
        }
        Global.arrayCurrentPopulation[p].mutate=true;

        for (int i=3;i<Global.seqLength-3;i++){
            c=_operator.mutateConformation(Global.arrayCurrentPopulation[p], i);
            if (c!=null) {
                neighbor.add(new Conformation(c));
            }
          }
        Collections.sort(neighbor);
        if (!Common.isExist(neighbor.get(0),pop,_pctAllowautoTryCount)) pop.add(neighbor.get(0));
    }
    Collections.sort(pop);

    int i=0;
    for (Conformation el: pop){
        Global.arrayCurrentPopulation[i++]=new Conformation(el);
        if (i>=Global.populationSize) break;
    }
}

Conformation reverseConformation(Conformation _c){
    AminoAcid[] aa=new AminoAcid[seqLen];
    int k=0;
    for (int i=seqLen-1;i>-1;i--){
        aa[k++]=new AminoAcid(_c.aminoAcid[i]);
    }
    RelativeZone rz=new RelativeZone();
    String sequence=rz.encode(aa);
    int[] fitnes=new int[2];
    Common.calculateFitness(aa,fitnes);
    Conformation c = new Conformation(aa, fitnes, sequence);
    return c;
}


public void initPopulationGenerate() throws Exception{
    Conformation cx=new Conformation();
    //ArrayList<Conformation> temp=new ArrayList<Conformation>();

    //Conformation[] temp=new Conformation[Global.reservedSize];
    boolean isExist=false;
    int size=0;
    Global.alReservedSol.clear();
    size=Global.alReservedSol.size();
    while (size<Global.reservedSize){
        Conformation c =cx.randConformation();
        for(Conformation el:Global.alReservedSol){
            if (el.encodedString.equalsIgnoreCase(c.encodedString)){
               isExist=true;
               break;
            }
        }
        if (!isExist){
            isExist=false;

            Global.alReservedSol.add(c);
            size=Global.alReservedSol.size();
        }
    }
    //Common.sortByFitness(Global.alReservedSol);
    Collections.sort(Global.alReservedSol);

    for (int i = 0; i < Global.populationSize; i++) {
        if(i%2==1){
                Conformation c=reverseConformation(Global.alReservedSol.get(i));
                Global.arrayInitialPopulation[i]=c;
                Global.arrayCurrentPopulation[i]=c;
                continue;
        }
        Global.arrayInitialPopulation[i]=Global.alReservedSol.get(i);
        Global.arrayCurrentPopulation[i]=Global.alReservedSol.get(i);

    }
}

}
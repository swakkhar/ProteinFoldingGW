/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ga;

import java.util.Random;
import java.util.ArrayList;
import java.io.*;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.util.Scanner;
/**
 *
 * @author mar
 */
public class Main {
    static String tag="";
    static int jobs=0;
    static int start=0;
    static int end=0;
    static String sOutput="";

    static String outputFile;
    static String baseDir;
    static boolean isHESF=false;
    static int timeSec=1;
    static long rndSeed=0;
    static int tagIdx=0;
    static String strSeq="";
    static int run=0;

    //this param need to be tested
    static int _pctAllowSimilar=100;

    private static void initGlobalParams(){
        Global.arrayInitialPopulation=new Conformation[Global.populationSize];
        Global.arrayCurrentPopulation=new Conformation[Global.populationSize];
    }
    private static void printCords(String cordFile) throws FileNotFoundException, IOException{
        FileOutputStream fCord=new FileOutputStream (cordFile);
        new PrintStream(fCord).println ("--Structure: "+Global.bestFitConf.encodedString);
        new PrintStream(fCord).println ("--Energy: "+Global.bestFitness+", No. of iterations: "+Global.maxTry +" ******");
        new PrintStream(fCord).println ("-- X  "+" "+"Y  "+" "+"Z  ");
        new PrintStream(fCord).println ("----------------------------------");
        for (AminoAcid el:Global.bestFitConf.aminoAcid){
            new PrintStream(fCord).println ((200+el.x)+" "+(200+el.y)+" "+(200+el.z));
        }
        new PrintStream(fCord).println ("----------------------------------");
        fCord.close();
    }

    public static void main(String[] args) throws Exception{
        int noofrun=55;
        tagIdx=0;
        if(args.length>2){
            tagIdx=Integer.parseInt(args[2])-1;
        }
        if(args.length>1){
            noofrun=Integer.parseInt(args[1]);
        }
        noofrun+=tagIdx;

        while(tagIdx<noofrun){
            String strSeq="";
            String cordFile="";
            String logFile="";
            String oplogFile="";
            Calendar currentDate = Calendar.getInstance();
            SimpleDateFormat formatter= new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            SimpleDateFormat formatter2= new SimpleDateFormat("yyyyMMdd-HHmmssa");
            String dateNow = formatter.format(currentDate.getTime());
            String fid=formatter2.format(currentDate.getTime());

            baseDir="output/";
            run=101;
            tagIdx++;

            int benchmarkindex=16;
            if(args.length>0){
                benchmarkindex=Integer.parseInt(args[0]);
            }
            if(benchmarkindex<0||benchmarkindex>16){
                System.out.println("No such benchmark");
                return;
            }
            

            strSeq=Global.benchmarkprimary[benchmarkindex][1].toUpperCase();
            tag=Global.benchmarkprimary[benchmarkindex][0];//"4PTI";

            System.out.println(Global.benchmarkprimary[benchmarkindex][1]);
            System.out.println(Global.benchmarkprimary[benchmarkindex][0]);


            Global.populationSize=50;

            Global.maxTry=Global.maxTry=(int)(50.0/Global.populationSize*1000*Global.iterationparam[benchmarkindex]);            

            if(args.length>3){
                Global.populationSize=Integer.parseInt(args[3]);
                Global.maxTry=(int)(50.0/Global.populationSize*1000*Global.iterationparam[benchmarkindex]);
            }

            Global.seqgreedlevel=Global.greedparam[benchmarkindex];            
            

            Global.seqstagparam=(int)((50.0/Global.populationSize)*Global.stagparam[benchmarkindex]);
            Global.distswitchparam=(int)(Global.maxTry/50.0);
            

            Global.traceOn=false;
            Global.logging=true;

            Global.pid=tag;//+"_"+String.format("%03d",run);
            outputFile=baseDir+Global.pid+"_"+tagIdx+".txt";
            cordFile=baseDir+Global.pid+"_"+tagIdx+"c.txt";
            logFile=baseDir+Global.pid+"_"+tagIdx+"log.txt";
            oplogFile=baseDir+Global.pid+"_"+tagIdx+"logop.txt";

            initGlobalParams();
            //Global.traceOn=true;
            Global.fout = new FileOutputStream (outputFile);
            Global.foutlog = new FileOutputStream (logFile);
            Global.foutoplog = new FileOutputStream (oplogFile);

            sOutput="--********** Date & Time: "+ dateNow+", No. of iterations: "+Global.maxTry;
            new PrintStream(Global.fout).println (sOutput);
            new PrintStream(Global.foutlog).println (sOutput);
            sOutput="--Sequence: "+strSeq;
            new PrintStream(Global.fout).println (sOutput);
            new PrintStream(Global.foutlog).println (sOutput);

           Random rnd=new Random();
           Global.rndSeed=Math.abs(rnd.nextInt());
           Global.ranDOM=new Random(Global.rndSeed);

           char a,b;
           Scanner energyRead = new Scanner(new BufferedReader(new FileReader("table.txt")));
           for(int i=0; i<400;i++){
               a=energyRead.next().toUpperCase().charAt(0);
               b=energyRead.next().toUpperCase().charAt(0);
               Global.energyfunction[a][b]= Integer.parseInt(energyRead.next());
               //System.out.println(a+" "+b+" "+Global.energyfunction[a][b]);
           }
           energyRead.close();

           GAPlus ga= new GAPlus(strSeq);

           //depends on ga constructor
           new EnergyDistribution();

           ga.initPopulationGenerate();

           sOutput="--rndSeed: "+Global.rndSeed+" popSize: "+Global.populationSize+", seqLen: "+Global.seqLength+", pctSimilar: "+_pctAllowSimilar+" *********************";
           new PrintStream(Global.foutlog).println (sOutput);
           sOutput="--***********************************************************************************************";
           new PrintStream(Global.foutlog).println (sOutput);

           ga.startSearch(_pctAllowSimilar);

           Global.fout.close();
           Global.foutlog.close();
           Global.foutoplog.close();

           //print coordinates
           printCords(cordFile);
        }
    }
}
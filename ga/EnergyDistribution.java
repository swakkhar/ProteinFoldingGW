/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author Ahammed Ullah
 */
public class EnergyDistribution {

    public EnergyDistribution() {

        ArrayList<EnergyMapping> acidenergymap=new ArrayList<EnergyMapping>();
        ArrayList<EnergyMapping> acidenergymap1=new ArrayList<EnergyMapping>();

        for(int i=0; i<Global.seqLength-1; i++){
            for(int j=i+2; j<Global.seqLength; j++){
                acidenergymap.add(new EnergyMapping(new AcidPair(i, j),
                        Global.energyfunction[Global.aaseq.charAt(i)][Global.aaseq.charAt(j)]));

                //real distribution assignment
                Global.energydistributions[0][i][j]=Global.energyfunction[Global.aaseq.charAt(i)][Global.aaseq.charAt(j)];
            }
        }
        Collections.sort(acidenergymap);
        //System.out.println(acidenergymap.size());

        //System.out.println("Original energy distribution");
        //printDistSorted(acidenergymap);

        acidenergymap1=H1MMJ(acidenergymap);
        assignSimulatedDist(acidenergymap1);
        //print
        //System.out.println("Simulated energy distribution");
        printDistSorted(acidenergymap1,acidenergymap);
    }
    private void assignSimulatedDist(ArrayList<EnergyMapping> acidenergysim){
        for(int i=0; i<acidenergysim.size(); i++){
            Global.energydistributions[1][acidenergysim.get(i).acidab.acida]
                    [acidenergysim.get(i).acidab.acidb]=acidenergysim.get(i).energy;

            if(Global.traceOn){
                //System.out.print(acidenergysim.get(i).acidab.toString()+" ");
                //System.out.print(Global.energydistributions[1]
                  //      [acidenergysim.get(i).acidab.acida][acidenergysim.get(i).acidab.acidb]+"\n");
            }

        }
    }

    private ArrayList<EnergyMapping> H1MMJ(ArrayList<EnergyMapping> acidenergymap){
        ArrayList<EnergyMapping> acidenergymap1=new ArrayList<EnergyMapping>();

        int saveindex=0;
        for(int i=0; i<acidenergymap.size(); i++){
            if(acidenergymap.get(i).energy>0){
                saveindex=i;
                break;
            }
        }

        AcidPair apair;

        int incrfactcmn=Global.seqgreedlevel;
        int weight=1;
        int incrementfact=incrfactcmn;
        int increment=incrementfact;

        int greeddivindex=0;
        int weightintoenergy;

        int ncontactsame=0;
        for(int i=saveindex; i<acidenergymap.size(); i++){
            ncontactsame++;
            while(Math.abs(acidenergymap.get(i).energy)>Global.greeddivision[greeddivindex]){
                //if(Global.traceOn)
                  //  System.out.println(acidenergymap.get(i-1).acidab.toString()+" "+weight+" :"+ncontactsame);
                weight+=increment;
                increment+=incrementfact;
                ncontactsame=0;

                greeddivindex++;
            }

            apair=acidenergymap.get(i).acidab;
            weightintoenergy=weight*acidenergymap.get(i).energy;//new
            acidenergymap1.add(new EnergyMapping(apair,weightintoenergy));//new
        }
        ncontactsame++;
        //if(Global.traceOn)
           // System.out.println(acidenergymap.get(acidenergymap.size()-1).acidab.toString()+" "+weight+" :"+ncontactsame);


        weight=1;
        incrementfact=incrfactcmn;
        increment=incrementfact;

        greeddivindex=0;
        ncontactsame=0;

        for(int i=saveindex-1; i>=0; i--){
            ncontactsame++;
            while(Math.abs(acidenergymap.get(i).energy)>Global.greeddivision[greeddivindex]){
                //if(Global.traceOn)
                  //  System.out.println(acidenergymap.get(i+1).acidab.toString()+" "+weight+" :"+ncontactsame);
                weight+=increment;
                increment+=incrementfact;
                ncontactsame=0;
                greeddivindex++;
            }

            apair=acidenergymap.get(i).acidab;
            weightintoenergy=weight*acidenergymap.get(i).energy;//new
            acidenergymap1.add(new EnergyMapping(apair,weightintoenergy));//new
        }
        ncontactsame++;
        //if(Global.traceOn)
        //    System.out.println(acidenergymap.get(0).acidab.toString()+" "+weight+" :"+ncontactsame);

        Collections.sort(acidenergymap1);
        return acidenergymap1;
    }

    public void printDistSorted(ArrayList<EnergyMapping> acidenergymap,ArrayList<EnergyMapping> real){
        int counter=0;
        for(int i=0; i<acidenergymap.size()-1; i++){
            counter++;
            if(acidenergymap.get(i).energy!=acidenergymap.get(i+1).energy){
                //System.out.print(acidenergymap.get(i).acidab.toString()+" "+acidenergymap.get(i).energy);
                System.out.print(acidenergymap.get(i).energy/real.get(i).energy+" "+acidenergymap.get(i).energy+" "+real.get(i).energy);
                System.out.println(" "+counter);
                counter=0;
            }
        }
        counter++;
        //System.out.print(acidenergymap.get(acidenergymap.size()-1).acidab.toString()+" "+acidenergymap.get(acidenergymap.size()-1).energy);
        System.out.print(acidenergymap.get(acidenergymap.size()-1).energy/real.get(acidenergymap.size()-1).energy+" "+acidenergymap.get(acidenergymap.size()-1).energy+" "+real.get(acidenergymap.size()-1).energy);
        System.out.println(" "+counter);
    }


    class AcidPair{
        int acida;
        int acidb;

        public AcidPair(int acida, int acidb) {
            this.acida = acida;
            this.acidb = acidb;
        }

        public boolean equals(Object ab){
            AcidPair p=(AcidPair)ab;
            return (acida == p.acida) && (acidb == p.acidb);
        }
        public int hashCode(){
            return (acida+" "+acidb).hashCode();
        }
        public String toString(){
            return "(" + acida + ", " + acidb + ")";
        }
    }

    class EnergyMapping implements Comparable{
        public AcidPair acidab;
        public int energy;

        public EnergyMapping(AcidPair acidab, int energy) {
            this.acidab =new AcidPair(acidab.acida, acidab.acidb);
            this.energy = energy;
        }

        public int compareTo(Object ob){
            if(this.energy == ((EnergyMapping) ob).energy)
                return 0;
            else if(this.energy > ((EnergyMapping) ob).energy)
                return 1;
            else
                return -1;
        }
    }
}
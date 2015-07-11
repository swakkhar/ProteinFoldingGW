/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

/**
 *
 * @author mar
 */
import java.util.ArrayList;
//import java.util.Random;
//import java.util.Enumeration;
import java.io.*;

import java.io.BufferedReader;
import java.io.IOException;
public class Common {
    Operators _operator;
    public Common() {
        _operator=new Operators();
    }
        
public static boolean isExist(Conformation _c, ArrayList<Conformation> _list,int _pctSimlilar){
    int _similar=100;
    if (_pctSimlilar>=99){
        for(Conformation el:_list){
            if (el.encodedString.equalsIgnoreCase(_c.encodedString)){
                Global.discardedTwins++;
                return true;
            }
        }
    }else{
        for(Conformation el:_list){
            _similar=getSimilarity(el.encodedString, _c.encodedString);
            if (_similar>=_pctSimlilar){
               return true;
            }
        }
    }
    return false;
}
public static int getRandomNext(int start,int end){     
     int next=(Math.abs(Global.ranDOM.nextInt())%(end-start))+start;
     return next;
 }

public static int checkValidity(AminoAcid[] aa) {
        int result = 0; //0 is valid; 1 is overlapped; 2 is distance is not equal to sqrt(2)
        int dx = 0, dy = 0, dz = 0;
        if (aa[Global.seqLength-1]==null)
            Global.SomethingWrong=true;
        for (int i = 1; i < Global.seqLength; i++) {
            dx = aa[i].x - aa[i - 1].x;
            dy = aa[i].y - aa[i - 1].y;
            dz = aa[i].z - aa[i - 1].z;
            if (Math.abs(dx)>1 || Math.abs(dy)>1 || Math.abs(dz)>1) {
                result = 2;
                return result;
            }
            if (Math.abs(dx) + Math.abs(dy) + Math.abs(dz) != 2) {
                result = 2;
                return result;
            }
        }
        //self-avoiding-walk
        int x = 0, y = 0, z = 0, x1 = 0, y1 = 0, z1 = 0;
        for (int i = 0; i < Global.seqLength; i++) {
            x = aa[i].x;
            y = aa[i].y;
            z = aa[i].z;
            for (int j =i+1; j < Global.seqLength; j++) {
//                for (int j =0; j < Global.seqLength; j++) {
//                if (i!=j){
                x1 = aa[j].x;
                y1 = aa[j].y;
                z1 = aa[j].z;
                if (x == x1 && y == y1 && z == z1) {
//                    System.out.println("{"+x+","+y+","+z+"} , {"+x1+","+y1+","+z1+"}");
                    result = 1;
                    return result;
                }
//                }
            }
        }

        return result;
    }

public static void calculateFitness(AminoAcid[] aa,int fitness[]) {
    fitness[0]=0;
    fitness[1]=0;
    
    int dx,dy,dz,d2;
    for(int i=0; i<aa.length-1; i++){
        for(int j=i+2; j<aa.length; j++){
            dx=aa[i].x-aa[j].x;
            dy=aa[i].y-aa[j].y;
            dz=aa[i].z-aa[j].z;
            d2=dx*dx+dy*dy+dz*dz-2;
            if (d2==0) {
                //energy*1000 form                
                fitness[0]+=Global.energydistributions[0][i][j];
                //simulated fitness form(search driver)
                fitness[1]+=Global.energydistributions[Global.energydistselected][i][j];
            }
        }
    }    
}
public static int getSimilarity(String S1,String S2){
    int len=S1.length();
    String s1;
    String s2;
    int similar=0,pct=100;
    if (S1.equalsIgnoreCase(S2)){
        return pct;
    }
    for(int i=0;i<len;i++){
        s1=S1.substring(i, i+1);
        s2=S2.substring(i, i+1);
        if(s1.equalsIgnoreCase(s2)){
            similar++;
        }
    }
    pct=(int)((1.0*similar+1.0)*100)/len;
    return pct;
}
public static ArrayList<Node3d> findCommonNeighbors(Node3d n1,Node3d n3){
    //n2 is the target point
    Node3d [] neigh_of_n1=new Node3d[12];
    Node3d [] neigh_of_n3=new Node3d[12];

    neigh_of_n1=Common.generateTN(n1);
    neigh_of_n3=Common.generateTN(n3);

    ArrayList<Node3d> commNeighbors=new ArrayList<Node3d>();

    for(Node3d e1:neigh_of_n1){
        for(Node3d e3:neigh_of_n3){
            if(isEqual(e1,e3)){
                commNeighbors.add(new Node3d(e1));
            }
        }
    }
    if (commNeighbors.size()<1){
        Global.SomethingWrong=false;
    }
    //if (Global.traceOn) System.out.println(commNeighbors.size());
    return commNeighbors;
}
    public static boolean isEqual(Node3d node1,Node3d node2){
        if (node1.x!=node2.x || node1.y!=node2.y || node1.z!=node2.z)
            return false;
        else
        return true;
    }

    public static boolean isTopologicalNeighbor(Node3d node1, Node3d node2) {
        if (getDistanceSqr(node1, node2)==2)
            return true;
        else
            return false;
    }

    public static boolean isFree(Node3d targetNode,ArrayList<Node3d> _occupiedNode) {
        for (Node3d el : _occupiedNode) {
            if (el.x == targetNode.x && el.y == targetNode.y && el.z == targetNode.z) {
                return false;
            }
        }
        return true;
    }
     public static boolean isFree(Node3d targetNode,Conformation _c) {
        for (AminoAcid el : _c.aminoAcid) {
            if (el.x == targetNode.x && el.y == targetNode.y && el.z == targetNode.z) {
                return false;
            }
        }
        return true;
    }

    public static Node3d[] generateTN(Node3d n3d) {
        int x = n3d.x;
        int y = n3d.y;
        int z = n3d.z;
        Node3d[] neighbor=new Node3d[12];
        neighbor[0] = new Node3d(x + 1, y + 1, z + 0);
       neighbor[1] = new Node3d(x + 1, y - 1, z + 0);
       neighbor[2] = new Node3d(x - 1, y - 1, z + 0);
       neighbor[3] = new Node3d(x - 1, y + 1, z + 0);
        //yz-plane @x=0 and denoted the vectors as e,f,g, and h
       neighbor[4] = new Node3d(x + 0, y + 1, z + 1);
       neighbor[5] = new Node3d(x + 0, y - 1, z + 1);
       neighbor[6] = new Node3d(x + 0, y - 1, z - 1);
       neighbor[7] = new Node3d(x + 0, y + 1, z - 1);
        //zx-plane @y=0 and denoted the vectors as i,j,k, and l
       neighbor[8] = new Node3d(x + 1, y + 0, z + 1);
       neighbor[9] = new Node3d(x + 1, y + 0, z - 1);
       neighbor[10] = new Node3d(x - 1, y + 0, z - 1);
       neighbor[11] = new Node3d(x - 1, y + 0, z + 1);
       return neighbor;
    }
    public static char getDirection(Node3d node1, Node3d node2) {
        char str = 0;
        int dx, dy, dz;

        dx = node2.x - node1.x;
        dy = node2.y - node1.y;
        dz = node2.z - node1.z;

        for (int i = 0; i < 12; i++) {
            if (dx == Global.vector3Dfcc[i].x & dy == Global.vector3Dfcc[i].y & dz == Global.vector3Dfcc[i].z) {
                int ascii = i + 97;
                str = (char) ascii;
                break;
            }
        }

        return str;
    }
public static ArrayList<Node3d> generateFreeTN(Conformation c, Node3d n3d){
    int x = n3d.x;
    int y = n3d.y;
    int z = n3d.z;
    Node3d[] neighbor=new Node3d[12];
    ArrayList<Node3d> tempList=new ArrayList<Node3d>();

    neighbor[0] = new Node3d(x + 1, y + 1, z + 0);
    neighbor[1] = new Node3d(x + 1, y - 1, z + 0);
    neighbor[2] = new Node3d(x - 1, y - 1, z + 0);
    neighbor[3] = new Node3d(x - 1, y + 1, z + 0);
    //yz-plane @x=0 and denoted the vectors as e,f,g, and h
    neighbor[4] = new Node3d(x + 0, y + 1, z + 1);
    neighbor[5] = new Node3d(x + 0, y - 1, z + 1);
    neighbor[6] = new Node3d(x + 0, y - 1, z - 1);
    neighbor[7] = new Node3d(x + 0, y + 1, z - 1);
    //zx-plane @y=0 and denoted the vectors as i,j,k, and l
    neighbor[8] = new Node3d(x + 1, y + 0, z + 1);
    neighbor[9] = new Node3d(x + 1, y + 0, z - 1);
    neighbor[10] = new Node3d(x - 1, y + 0, z - 1);
    neighbor[11] = new Node3d(x - 1, y + 0, z + 1);
    for (Node3d n:neighbor){
        if (isFree(n, c)){
            tempList.add(n);
        }
    }
    return tempList;
}
   public static char newDirection(char chr) {
        char ch = chr;
        ////Random r = new Random();
        int i = Common.getRandomNext(0,12);
        int ascii = i + 97;
        if ((char) ascii != chr) {
            ch = (char) ascii;
        } else {
            newDirection(ch);
        }
        return ch;
    }

    public static double getDistance(Node3d point1, Node3d point2){
        double d=0;
        int dx=point2.x-point1.x;
        int dy=point2.y-point1.y;
        int dz=point2.z-point1.z;
        d=Math.sqrt(dx*dx+dy*dy+dz*dz);
        return d;
    }
    public static int getDistanceSqr(Node3d point1, Node3d point2){
        int d=0;
        int dx=point2.x-point1.x;
        int dy=point2.y-point1.y;
        int dz=point2.z-point1.z;
        d=dx*dx+dy*dy+dz*dz;
        return d;
    }
    public static void  sortByFitness(ArrayList<Conformation> popList) throws Exception{
        int size=popList.size();
        Conformation temp=new Conformation();
        
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                if (popList.get(j).fitness[1]<popList.get(i).fitness[1]) {
                    temp=popList.get(i);
                    popList.set(i, popList.get(j));
                    popList.set(j,temp);
                }
            }
        }
    }
    public static void  sortByFitness(Conformation[] popList) throws Exception{
        int size=popList.length;
        Conformation temp=new Conformation();

        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                if (popList[j].fitness[1]<popList[i].fitness[1]) {
                    temp=popList[i];
                    popList[i]= popList[j];
                    popList[j]=temp;
                }
            }
        }
    }
}
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
import java.util.Collections;
//import sun.awt.SunToolkit.OperationTimedOut;

public class Operators{
private ArrayList<Node3d> findLxCx(Node3d ni, Node3d ni1,Conformation _c){
        //Lx is the new point where to pull
        //Cx is the common neighbor of Lx and point to be pulled (ni)
        //ni1 is either n(i-1) or n(i+1)

        Node3d [] neigh_of_ni1=new Node3d[12];
        Node3d [] neigh_of_Lx=new Node3d[12];
        Node3d Lx=new Node3d();
        Node3d Cx=new Node3d();

        ArrayList<Node3d> possibleLx=new ArrayList<Node3d>();
        ArrayList<Node3d> LxCx=new ArrayList<Node3d>();
        int x,y,z;
        neigh_of_ni1=Common.generateTN(ni1);
        for(int i=0;i<12;i++){
            Lx=new Node3d(neigh_of_ni1[i]);
            if(Common.isFree(Lx,_c)){
               possibleLx.add(Lx);
            }
        }

        for (Node3d lx:possibleLx){
            neigh_of_Lx=Common.generateTN(lx);
            for (Node3d cx:neigh_of_Lx){
                if(Common.isFree(cx,_c)){
                    if(Common.isTopologicalNeighbor(cx, ni)){
                        LxCx.add(lx);
                        LxCx.add(cx);
                        return LxCx;
                    }
                }
            }

        }
        return LxCx;
    }
public Conformation mutateConformation(Conformation _c, int _pos) throws Exception{
        Conformation c;// = new Conformation(Ci.aminoAcid, Ci.fitness, Ci.encodedString);

        String str = "";
        RelativeZone rz=new RelativeZone();
        str = _c.encodedString;
        char[] ch = str.toCharArray();
        char chr = Common.newDirection(ch[_pos]);
        if ((int) chr == 0) {
            c=new Conformation(_c);
        } else {
            ch[_pos] = chr;
            str="";
            for (char el:ch){
                str+=Character.toString(el);
            }
            c=rz.decode(str);
            if (c==null)
                return null;
        }
        c.aminoAcid[_pos].operation[OperatorType.CLASSIC_MUTATION]=Global.generationCount;
        Global.solExplored+=1;
        return c;
    }
public void crossOver(int idx1,int idx2,int pos, Conformation[] outC) throws Exception{
    //Conformation c = new Conformation();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    RelativeZone rz=new RelativeZone();
    int size=Global.seqLength;
    boolean succeed=false;
    String s1,s2,s1p1,s1p2,s2p1,s2p2;
    int si,si2;

    Conformation c1=new Conformation(Global.arrayCurrentPopulation[idx1]);
    Conformation c2=new Conformation(Global.arrayCurrentPopulation[idx2]);
//
//    Operation[] op1=new Operation[size];
//    Operation[] op2=new Operation[size];
//    conf1.copyOperation(op1);
//    conf2.copyOperation(op2);

    s1=c1.encodedString;
    si=s1.length();
    s2=c2.encodedString;
    si2=s2.length();
    if (si==si2 & si<size-2)
        Global.SomethingWrong=true;

    s1p1=s1.substring(0, pos);
    s1p2=s1.substring(pos);
    s2p1=s2.substring(0, pos);
    s2p2=s2.substring(pos);

    Global.xoverAttempts+=1;
    //System.out.println(Global.xoverAttempts);
    s1=s1p1+s2p2;
    s2=s2p1+s1p2;
//    try{

        ArrayList<Conformation> tempList=new ArrayList<Conformation>();

//        tempList.add(c1);
//        tempList.add(c2);

        Conformation c1x=rz.decode(s1);

        if (c1x!=null && Common.checkValidity(c1x.aminoAcid)==0){
            Global.solExplored+=1;
            int i=0;
            while(i<pos){
                //c1x.operation[i]=new Operation(op1[i]);
                c1x.aminoAcid[i].operation=c1.aminoAcid[i].operation;
                i++;
            }
            i=pos;
            int []op={-1,-1,-1,-1,-1};
            op[OperatorType.CROSS_OVER]=Global.generationCount;
            c1x.aminoAcid[i].operation=op;
            i++;
            while(i<size){
                //c1x.operation[i]=new Operation(op2[i]);
                c1x.aminoAcid[i].operation=c2.aminoAcid[i].operation;
                i++;
            }
            succeed=true;
            tempList.add(c1x);
         }
        Conformation c2x=rz.decode(s2);
        if (c2x!=null && Common.checkValidity(c2x.aminoAcid)==0){
           Global.solExplored+=1;
           int i=0;
            while(i<pos){
                //c2x.operation[i]=new Operation(op2[i]);
                c2x.aminoAcid[i].operation=c2.aminoAcid[i].operation;
                i++;
            }
            i=pos;
            int []op={-1,-1,-1,-1,-1};
            op[OperatorType.CROSS_OVER]=Global.generationCount;
            c2x.aminoAcid[i].operation=op;
            i++;
            while(i<size){
                //c2x.operation[i]=new Operation(op1[i]);
                c2x.aminoAcid[i].operation=c1.aminoAcid[i].operation;
                i++;
            }

            succeed=true;
           tempList.add(c2x);
         }
        if (succeed) {
            Global.xoverSucceeds+=1;
             Collections.sort(tempList);
             outC[0]=new Conformation(tempList.get(0));
            if (tempList.size()>1) outC[1]=new Conformation(tempList.get(1));
        }
    else
     outC=null;
    }
public Conformation pullMove(Conformation _c,int _posx) throws Exception{
    Node3d n1,n,n2,_n,temp;//n1=_pos-1,n=_pos,n2=_pos+1,_n=new pos
    int x,y,z,dir, size,_pos;
    boolean isTN1=false,isTN2=false,isTN12=false;
    ArrayList<Node3d> freeNode=new ArrayList<Node3d>();
    ArrayList<Node3d> LxCx=new ArrayList<Node3d>();
    Node3d Cx=new Node3d();
    Node3d Lx=new Node3d();

    Global.pullAttempts+=1;

    size=Global.seqLength;
    AminoAcid aa[]=new AminoAcid[size];

    if (_posx==0)
        _pos=1;
    else if (_posx==size-1)
        _pos=_posx-1;
    else
        _pos=_posx;

    n1=new Node3d(_c.aminoAcid[_pos-1].x,_c.aminoAcid[_pos-1].y,_c.aminoAcid[_pos-1].z);
    n= new Node3d(_c.aminoAcid[_pos].x,_c.aminoAcid[_pos].y,_c.aminoAcid[_pos].z);
    n2=new Node3d(_c.aminoAcid[_pos+1].x,_c.aminoAcid[_pos+1].y,_c.aminoAcid[_pos+1].z);
    _n=new Node3d();

    if(Common.getRandomNext(0,2) == 0)
        isTN2=true;
     else
        isTN1 = true;

    if(isTN1){
        freeNode=getFreeNode(n1, _c);
        if (freeNode.size()>0){
            for(Node3d el:freeNode){
                if(Common.isTopologicalNeighbor(el, n2)){
                    isTN12=true;_n=new Node3d(el);
                    break;
                }
            }
            if (!isTN12){
                 LxCx=findLxCx(n, n1, _c);
            }
        }
    }
    else if(isTN2){
        freeNode=getFreeNode(n2, _c);
        if (freeNode.size()>0){
            for(Node3d el:freeNode){
                if(Common.isTopologicalNeighbor(el, n1)){
                    isTN12=true;_n=new Node3d(el);
                    break;
                }
            }
            if (!isTN12){
                LxCx=findLxCx(n, n2, _c);
            }
        }
     }
    if (!isTN12){
        if(LxCx.size()>1){
            Lx=LxCx.get(0);
            Cx=LxCx.get(1);
         }
        else{
            _c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
         return null;
         }
     }
    if(isTN12){
        _c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
         return null;
    }else if(isTN1){
        for(int i=0;i<_pos;i++) aa[i]=new AminoAcid(_c.aminoAcid[i]);
        aa[_pos]=new AminoAcid(Lx.x,Lx.y,Lx.z,_c.aminoAcid[_pos].operation);//Lx
        aa[_pos+1]=new AminoAcid(Cx.x,Cx.y,Cx.z,_c.aminoAcid[_pos].operation);//Cx
        //<pulling here>
        if (_pos<size-2){
            for( int count=_pos+2; count<size;count++){
                Node3d ni=new Node3d(_c.aminoAcid[count].x,_c.aminoAcid[count].y,_c.aminoAcid[count].z);
                Node3d ni1=new Node3d(aa[count-1].x,aa[count-1].y,aa[count-1].z);
                if(Common.isTopologicalNeighbor(ni, ni1)){
                    for (int j=count;j<size;j++){
                        aa[j]=new AminoAcid(_c.aminoAcid[j].x,_c.aminoAcid[j].y,_c.aminoAcid[j].z,_c.aminoAcid[j].operation);
                    }
                    break;
                }
                aa[count]=new AminoAcid(_c.aminoAcid[count-2].x,_c.aminoAcid[count-2].y,_c.aminoAcid[count-2].z,_c.aminoAcid[count-2].operation);
            }
        }
    }else if(isTN2){
        for(int i=size-1;i>_pos;i--) {
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
        aa[_pos]=new AminoAcid(Lx.x,Lx.y,Lx.z,_c.aminoAcid[_pos].operation);//Lx
        aa[_pos-1]=new AminoAcid(Cx.x,Cx.y,Cx.z,_c.aminoAcid[_pos].operation);//Cx
        //<pulling here>
        if (_pos>1){
            for( int count=_pos-2; count>-1;count--){
                Node3d ni=new Node3d(_c.aminoAcid[count].x,_c.aminoAcid[count].y,_c.aminoAcid[count].z);
                Node3d ni1=new Node3d(aa[count+1].x,aa[count+1].y,aa[count+1].z);
                if(Common.isTopologicalNeighbor(ni, ni1)){
                    for (int j=count;j>-1;j--){
                        aa[j]=new AminoAcid(_c.aminoAcid[j].x,_c.aminoAcid[j].y,_c.aminoAcid[j].z,_c.aminoAcid[j].operation);
                    }
                    break;
                }
                aa[count]=new AminoAcid(_c.aminoAcid[count+2].x,_c.aminoAcid[count+2].y,_c.aminoAcid[count+2].z,_c.aminoAcid[count+2].operation);
            }
        }
    }
    for(AminoAcid a:aa){
        if(a==null){
             _c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
            return null;
        }
            //Global.SomethingWrong=true;
    }

    if (Common.checkValidity(aa)==0){
        Global.solExplored+=1;
        RelativeZone rz=new RelativeZone();
        String str=rz.encode(aa);
        int fitness[]={0,0};Common.calculateFitness(aa,fitness);
        Conformation c= new Conformation(aa, fitness, str);

        c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
         return c;
     }
     else
    {
         _c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
         return null;
     }
 }
public Conformation pullMoveSinglePoint(Conformation _c,int _posx) throws Exception{
    Node3d curr_n,prev_n,next_n;//n1=_pos-1,n=_pos,n2=_pos+1,_n=new pos
    int x,y,z, size,_pos;
    ArrayList<Node3d> freeNode=new ArrayList<Node3d>();
    ArrayList<Conformation> neighbors=new ArrayList<Conformation>();

    size=Global.seqLength;
    AminoAcid aa[]=new AminoAcid[size];

    if (_posx==0)
        _pos=1;
    else if (_posx==size-1)
        _pos=_posx-1;
    else
        _pos=_posx;

    //prev_n=new Node3d(_c.aminoAcid[_pos-1].x,_c.aminoAcid[_pos-1].y,_c.aminoAcid[_pos-1].z);
    curr_n= new Node3d(_c.aminoAcid[_pos].x,_c.aminoAcid[_pos].y,_c.aminoAcid[_pos].z);
    //next_n=new Node3d(_c.aminoAcid[_pos+1].x,_c.aminoAcid[_pos+1].y,_c.aminoAcid[_pos+1].z);

    freeNode=getFreeNode(curr_n, _c);
    neighbors.add(_c);

    if (freeNode.size()>0){
        for(Node3d el:freeNode){
            //pull towards 0
                for (int pos=size-1;pos>=_pos;pos--){
                    aa[pos]=new AminoAcid(_c.aminoAcid[pos]);
                }
                x=el.x;y=el.y;z=el.z;
                aa[_pos-1]=new AminoAcid(x,y,z,_c.aminoAcid[_pos-1].operation);
                for (int pos=_pos-2;pos>-1;pos--){
                    aa[pos]=new AminoAcid(_c.aminoAcid[pos+1]);
                }
                if (Common.checkValidity(aa)==0){

                    Global.solExplored+=1;
                    RelativeZone rz=new RelativeZone();
                    String str=rz.encode(aa);
                    int fitness[]={0,0};Common.calculateFitness(aa,fitness);
                    Conformation c= new Conformation(aa, fitness, str);

                    //c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
                    if(!Common.isExist(c, neighbors,100)) neighbors.add(c);
                 }

            //pull towards size
                for (int pos=0;pos<=_pos;pos++){
                    aa[pos]=new AminoAcid(_c.aminoAcid[pos]);
                }
                x=el.x;y=el.y;z=el.z;
                aa[_pos+1]=new AminoAcid(el.x,el.y,el.z,_c.aminoAcid[_pos+1].operation);
                for (int pos=_pos+2;pos<size;pos++){
                    aa[pos]=new AminoAcid(_c.aminoAcid[_pos-1]);
                }
                if (Common.checkValidity(aa)==0){
                    Global.solExplored+=1;
                    RelativeZone rz=new RelativeZone();
                    String str=rz.encode(aa);
                    int fitness[]={0,0};Common.calculateFitness(aa,fitness);
                    Conformation c= new Conformation(aa, fitness, str);

                    //c.aminoAcid[_posx].operation[OperatorType.PULL_MOVE]=Global.generationCount;
                    if(!Common.isExist(c, neighbors,100)) neighbors.add(c);
                 }
        }
    }
    else{
        return null;
    }
    Collections.sort(neighbors);
    return neighbors.get(0);
 }
public Conformation randomWalk(Conformation _c,int _posx) throws Exception{
    Node3d n1,n,n2,_n,temp;//n1=_pos-1,n=_pos,n2=_pos+1,_n=new pos
    int seqLen,_pos;
    boolean isTN1=false,isTN2=false,isTN12=false;
    ArrayList<Node3d> freeNode=new ArrayList<Node3d>();
    ArrayList<Node3d> LxCx=new ArrayList<Node3d>();
    Node3d Cx=new Node3d();
    Node3d Lx=new Node3d();
    Global.pullAttempts+=1;
    seqLen=Global.seqLength;

    AminoAcid aa[]=new AminoAcid[seqLen];

    if (_posx==0)
        _pos=1;
    else if (_posx==seqLen-1)
        _pos=_posx-1;
    else
        _pos=_posx;
    _n=new Node3d();

    n1=new Node3d(_c.aminoAcid[_pos-1].x,_c.aminoAcid[_pos-1].y,_c.aminoAcid[_pos-1].z);
    n= new Node3d(_c.aminoAcid[_pos].x,_c.aminoAcid[_pos].y,_c.aminoAcid[_pos].z);
    n2=new Node3d(_c.aminoAcid[_pos+1].x,_c.aminoAcid[_pos+1].y,_c.aminoAcid[_pos+1].z);


    if(Common.getRandomNext(0,2) == 0)
        isTN2=true;
     else
        isTN1 = true;

    if(isTN1){
        freeNode=getFreeNode(n1, _c);
        if (freeNode.size()>0){
            for(Node3d el:freeNode){
                if(Common.isTopologicalNeighbor(el, n2)){
                    isTN12=true;_n=new Node3d(el);
                    break;
                }
            }
            if (!isTN12){
                 LxCx=findLxCx(n, n1, _c);
            }
        }
    }
    else if(isTN2){
        freeNode=getFreeNode(n2, _c);
        if (freeNode.size()>0){
            for(Node3d el:freeNode){
                if(Common.isTopologicalNeighbor(el, n1)){
                    isTN12=true;_n=new Node3d(el);
                    break;
                }
            }
            if (!isTN12){
                LxCx=findLxCx(n, n2, _c);
            }
        }
     }
    if (!isTN12){
        if(LxCx.size()>1){
            Lx=LxCx.get(0);
            Cx=LxCx.get(1);
         }
        else{
         return null;
         }
     }
    if(isTN12){

        for(int i=0;i<_pos;i++){
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
        aa[_pos]=new AminoAcid(_n.x,_n.y,_n.z);
        for(int i=_pos+1;i<seqLen;i++) {
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
    }else if(isTN1){
        for(int i=0;i<_pos;i++) aa[i]=new AminoAcid(_c.aminoAcid[i]);
        aa[_pos]=new AminoAcid(Lx.x,Lx.y,Lx.z);//Lx
        aa[_pos+1]=new AminoAcid(Cx.x,Cx.y,Cx.z);//Cx
        //<pulling here>
        if (_pos<seqLen-2){
            for( int count=_pos+2; count<seqLen;count++){
                Node3d ni=new Node3d(_c.aminoAcid[count].x,_c.aminoAcid[count].y,_c.aminoAcid[count].z);
                Node3d ni1=new Node3d(aa[count-1].x,aa[count-1].y,aa[count-1].z);
                if(Common.isTopologicalNeighbor(ni, ni1)){
                    for (int j=count;j<seqLen;j++){
                        aa[j]=new AminoAcid(_c.aminoAcid[j].x,_c.aminoAcid[j].y,_c.aminoAcid[j].z);
                    }
                    break;
                }
                aa[count]=new AminoAcid(_c.aminoAcid[count-2].x,_c.aminoAcid[count-2].y,_c.aminoAcid[count-2].z);
            }
        }
    }else if(isTN2){
        for(int i=seqLen-1;i>_pos;i--) {
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
        aa[_pos]=new AminoAcid(Lx.x,Lx.y,Lx.z);//Lx
        aa[_pos-1]=new AminoAcid(Cx.x,Cx.y,Cx.z);//Cx
        //<pulling here>
        if (_pos>1){
            for( int count=_pos-2; count>-1;count--){
                Node3d ni=new Node3d(_c.aminoAcid[count].x,_c.aminoAcid[count].y,_c.aminoAcid[count].z);
                Node3d ni1=new Node3d(aa[count+1].x,aa[count+1].y,aa[count+1].z);
                if(Common.isTopologicalNeighbor(ni, ni1)){
                    for (int j=count;j>-1;j--){
                        aa[j]=new AminoAcid(_c.aminoAcid[j].x,_c.aminoAcid[j].y,_c.aminoAcid[j].z);
                    }
                    break;
                }
                aa[count]=new AminoAcid(_c.aminoAcid[count+2].x,_c.aminoAcid[count+2].y,_c.aminoAcid[count+2].z);
            }
        }
    }
    for(AminoAcid a:aa){
        if(a==null){
            return null;
        }
            //Global.SomethingWrong=true;
    }

    if (Common.checkValidity(aa)==0){
        Global.solExplored+=1;
        RelativeZone rz=new RelativeZone();
        String str=rz.encode(aa);
        int fitness[]={0,0};Common.calculateFitness(aa,fitness);
        Conformation c= new Conformation(aa, fitness, str);
        return c;
     }
     else
    {
         return null;
     }
 }
public Conformation cornerFlip(Conformation _c,int _posx) throws Exception{//latest

    Node3d n1,n,n2,_n,temp;//n1=_pos-1,n=_pos,n2=_pos+1,_n=new pos
    int x,y,z,dir, seqLen,_pos;
    boolean isTN12=false;
    ArrayList<Node3d> freeNode=new ArrayList<Node3d>();

   Global.flipAttempts+=1;
    seqLen=Global.seqLength;
    AminoAcid aa[]=new AminoAcid[seqLen];

    if (_posx==0)
        _pos=1;
    else if (_posx==seqLen-1)
        _pos=_posx-1;
    else
        _pos=_posx;

    n1=new Node3d(_c.aminoAcid[_pos-1].x,_c.aminoAcid[_pos-1].y,_c.aminoAcid[_pos-1].z);
    n= new Node3d(_c.aminoAcid[_pos].x,_c.aminoAcid[_pos].y,_c.aminoAcid[_pos].z);
    n2=new Node3d(_c.aminoAcid[_pos+1].x,_c.aminoAcid[_pos+1].y,_c.aminoAcid[_pos+1].z);
    _n=new Node3d();

    freeNode=getFreeNode(n1, _c);
    if (freeNode.size()>0){
        for(Node3d el:freeNode){
            if(Common.isTopologicalNeighbor(el, n2)){
                isTN12=true;
                _n=new Node3d(el);
                break;
            }
        }
    }

    if(isTN12){
        /*if (Global.arraySeq[_posx]=='0' && Common.getRandomNext(0,5)<1){
            _c.aminoAcid[_posx].operation[OperatorType.CORNER_FLIP]=Global.generationCount;
            return null;
        }*/
        for(int i=0;i<_pos;i++){
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
        aa[_pos]=new AminoAcid(_n.x,_n.y,_n.z,_c.aminoAcid[_pos].operation);
        for(int i=_pos+1;i<seqLen;i++) {
            aa[i]=new AminoAcid(_c.aminoAcid[i]);
        }
    } else
    {
         _c.aminoAcid[_posx].operation[OperatorType.CORNER_FLIP]=Global.generationCount;
         return null;
     }
    if (Common.checkValidity(aa)==0){
        Global.flipSucceeds+=1;
        RelativeZone rz=new RelativeZone();
        String str=rz.encode(aa);
        int fitness[]={0,0};Common.calculateFitness(aa,fitness);
        Conformation c= new Conformation(aa, fitness, str);
        c.aminoAcid[_posx].operation[OperatorType.CORNER_FLIP]=Global.generationCount;
         return c;
     }
     else
    {
        _c.aminoAcid[_posx].operation[OperatorType.CORNER_FLIP]=Global.generationCount;
         return null;
     }
 }

public Conformation crankShaftRotation(Conformation _c,int _loopsize) throws Exception{//latest
      Node3d n1,n2,n3,n4;
      int x,y,z,fit;
      int k1,k2,k3,k4;
      int xy,yz,zx,pivotplane,movingplane;

      ArrayList<Node3d> freeTNn1=new ArrayList<Node3d>();
      ArrayList<Node3d> freeTNn4=new ArrayList<Node3d>();
      ArrayList<Conformation> neighborSol=new ArrayList<Conformation>();
      ArrayList<Node3d> commTNn1n3=new ArrayList<Node3d>();
      fit=_c.fitness[0];

      if (_loopsize==4){
        for (int i=0;i<Global.seqLength-5;i++){
            k1=i;k2=i+1;k3=i+2;k4=i+3;

            x=_c.aminoAcid[k1].x; y=_c.aminoAcid[k1].y; z=_c.aminoAcid[k1].z; n1=new Node3d(x, y, z);
            x=_c.aminoAcid[k4].x; y=_c.aminoAcid[i+3].y; z=_c.aminoAcid[k4].z;n4=new Node3d(x, y, z);

            if (Common.isTopologicalNeighbor(n1, n4)){
                x=_c.aminoAcid[k3].x; y=_c.aminoAcid[k3].y; z=_c.aminoAcid[k3].z; n2=new Node3d(x, y, z);
                x=_c.aminoAcid[k4].x; y=_c.aminoAcid[k4].y; z=_c.aminoAcid[k4].z; n3=new Node3d(x, y, z);
                freeTNn1=Common.generateFreeTN(_c,n1);
                if (freeTNn1.size()<1) continue;
                for (Node3d n:freeTNn1){
                    commTNn1n3=Common.findCommonNeighbors(n, n3);
                    if(commTNn1n3.size()<1) continue;
                    for(Node3d _n:commTNn1n3){
                        if (Common.isFree(_n, _c)){
                            n2=new Node3d(n);
                            n3=new Node3d(_n);
                            _c.aminoAcid[k2]=new AminoAcid(n2.x, n2.y, n2.z, _c.aminoAcid[k2].operation);
                            _c.aminoAcid[k3]=new AminoAcid(n3.x, n3.y, n3.z, _c.aminoAcid[k3].operation);
                            if (Common.checkValidity(_c.aminoAcid)==0){
                                Global.solExplored+=1;
                                RelativeZone rz=new RelativeZone();
                                String str=rz.encode(_c.aminoAcid);
                                int fitness[]={0,0};Common.calculateFitness(_c.aminoAcid,fitness);
                                Conformation c= new Conformation(_c.aminoAcid, fitness, str);
                                //c.aminoAcid[_posx].operation[OperatorType.CORNER_FLIP]=Global.generationCount;
                                neighborSol.add(c);
                                //if (Global.consoleTrace)  System.out.println("crankShaft succeed");
                             }
                        }
                    }
                }
            }

        }
        if (neighborSol.size()>0){
            Collections.sort(neighborSol);
            return neighborSol.get(0);
        }
         else{
            return _c;
         }
    }else{
            return _c;
      }
 }
public Conformation tiltMove(Conformation _c,int _pos, int _next ) throws Exception{

        if (_next==_pos) return _c;
        int x,y,z, size;

        size=Global.seqLength;

        Node3d node=new Node3d();
        Node3d nextnode=new Node3d();
        Node3d cx=new Node3d();
        Node3d lx=new Node3d();
        ArrayList<Node3d> lxcx=new ArrayList<Node3d>();
        ArrayList<Node3d> lnodes=new ArrayList<Node3d>();
        ArrayList<Node3d> cnodes=new ArrayList<Node3d>();

        ArrayList<Conformation> cout=new ArrayList<Conformation>();

        AminoAcid aa[]=new AminoAcid[size];
//        Conformation c1=new Conformation();

        x=_c.aminoAcid[_pos].x; y=_c.aminoAcid[_pos].y;  z=_c.aminoAcid[_pos].z;
        node=new Node3d(x,y,z);
        x=_c.aminoAcid[_next].x; y=_c.aminoAcid[_next].y; z=_c.aminoAcid[_next].z;
        nextnode=new Node3d(x,y,z);

        lxcx.clear();
        lxcx=findLCTilt(node, nextnode,_c);

        if(lxcx==null) return _c;
        int opt=lxcx.size();
        if (opt>0){
//            lx=lxcx.get(0);
//            cx=lxcx.get(1);
            for(int i=0;i<opt;i++){
                if (i%2==0) lnodes.add(lxcx.get(i));
                if (i%2==1) cnodes.add(lxcx.get(i));
            }

        }
        else
        {
            return _c;
        }
  for(int k=0;k<lnodes.size();k++){
      lx=lnodes.get(k);
      cx=cnodes.get(k);
        if(_next<_pos){

            aa[_pos]=new AminoAcid(lx.x,lx.y,lx.z,_c.aminoAcid[_pos].operation);
            for( int count=_pos+1; count<size;count++){
                aa[count]=new AminoAcid(_c.aminoAcid[count-1]);
            }
            aa[_next]=new AminoAcid(cx.x,cx.y,cx.z,_c.aminoAcid[_next].operation);
            for( int count=_next-1; count>=0;count--){
                aa[count]=new AminoAcid(_c.aminoAcid[count+1]);
            }

        }else{
            aa[_pos]=new AminoAcid(lx.x,lx.y,lx.z,_c.aminoAcid[_pos].operation);

            for( int count=_pos-1; count>=0;count--){
                aa[count]=new AminoAcid(_c.aminoAcid[count+1]);
            }
            aa[_next]=new AminoAcid(cx.x,cx.y,cx.z,_c.aminoAcid[_next].operation);
            for( int count=_next+1; count<size;count++){
                aa[count]=new AminoAcid(_c.aminoAcid[count-1]);
            }

        }
         if (Common.checkValidity(aa)==0){
             Global.solExplored+=1;
             RelativeZone rz=new RelativeZone();
             String str=rz.encode(aa);
            int fitness[]={0,0};Common.calculateFitness(aa,fitness);
            Conformation c=new Conformation(aa, fitness, str);
            c.aminoAcid[_pos].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            c.aminoAcid[_next].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            cout.add(c);
         }

        }
        if(cout.size()<1)
         {
            _c.aminoAcid[_pos].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            _c.aminoAcid[_next].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            return _c;
         }
        Collections.sort(cout);
        return cout.get(0);
    }
    private ArrayList<Node3d> findLCTilt(Node3d pos, Node3d next,Conformation _c){
        Node3d [] pos_neigh=new Node3d[12];
        Node3d [] next_neigh=new Node3d[12];
        ArrayList<Node3d> lcnodes=new ArrayList<Node3d>();
        pos_neigh=Common.generateTN(pos);
        next_neigh=Common.generateTN(pos);
        for(Node3d n1:pos_neigh){
            if(!Common.isFree(n1, _c)) continue;
            for (Node3d n2:next_neigh){
                if(!Common.isFree(n2, _c)) continue;
                if(Common.isTopologicalNeighbor(n1, n2)){
                    lcnodes.add(n1);
                    lcnodes.add(n2);
//                        return lcnodes;
                }
            }
        }
        return lcnodes;
    }
public void crossOverMulti(int idx1,int idx2, Conformation[] outC) throws Exception{
    //Conformation c = new Conformation();
    ArrayList<Conformation> neighbor=new ArrayList<Conformation>();
    RelativeZone rz=new RelativeZone();
    int size=Global.seqLength;
    boolean succeed=false;
    String s1,s2,s1p1,s1p2,s1p3,s2p1,s2p2,s2p3;
    int si,si2;

    Conformation c1=new Conformation(Global.arrayCurrentPopulation[idx1]);
    Conformation c2=new Conformation(Global.arrayCurrentPopulation[idx2]);

    int start=Common.getRandomNext(3, Global.seqLength-3);
    int end=Common.getRandomNext(3, Global.seqLength-3);
    int temp;
    if (start>end) {
        temp=start;
        start=end;
        end=temp;
    }
    s1=c1.encodedString;
    si=s1.length();
    s2=c2.encodedString;
    si2=s2.length();
    if (si==si2 & si<size-2)
        Global.SomethingWrong=true;

    s1p1=s1.substring(0, start);
    s1p2=s1.substring(start,end);
    s1p3=s1.substring(end);

    s2p1=s2.substring(0, start);
    s2p2=s2.substring(start,end);
    s2p3=s2.substring(end);

    Global.xoverAttempts+=1;

    s1=s1p1+s2p2+s1p3;
    s2=s2p1+s1p2+s2p3;
//    try{

        ArrayList<Conformation> tempList=new ArrayList<Conformation>();

//        tempList.add(c1);
//        tempList.add(c2);

        Conformation c1x=rz.decode(s1);

        if (c1x!=null && Common.checkValidity(c1x.aminoAcid)==0){
            succeed=true;
            tempList.add(c1x);
         }
        Conformation c2x=rz.decode(s2);
        if (c2x!=null && Common.checkValidity(c2x.aminoAcid)==0){
           succeed=true;
           tempList.add(c2x);
         }
        if (succeed) {
            Global.xoverSucceeds+=1;
             Collections.sort(tempList);
             outC[0]=new Conformation(tempList.get(0));
            if (tempList.size()>1) outC[1]=new Conformation(tempList.get(1));
        } else
            outC=null;
    }
public Conformation tiltMove2(Conformation _c,int _pos, int _next ) throws Exception{

        if (_next==_pos) return _c;
        int x,y,z, size;

        size=Global.seqLength;

        Node3d node=new Node3d();
        Node3d nextnode=new Node3d();
        Node3d Cx=new Node3d();
        Node3d Lx=new Node3d();
        ArrayList<Node3d> LxCx=new ArrayList<Node3d>();

        AminoAcid aa[]=new AminoAcid[size];

        x=_c.aminoAcid[_pos].x; y=_c.aminoAcid[_pos].y;  z=_c.aminoAcid[_pos].z;
        node=new Node3d(x,y,z);
        x=_c.aminoAcid[_next].x; y=_c.aminoAcid[_next].y; z=_c.aminoAcid[_next].z;
        nextnode=new Node3d(x,y,z);

        LxCx.clear();
        LxCx=findLxCx(node, nextnode,_c);
        if (LxCx==null||LxCx.size()<2){
            //un successfull for this position
             _c.aminoAcid[_pos].operation[OperatorType.TILT_MOVE]= Global.generationCount;
             _c.aminoAcid[_next].operation[OperatorType.TILT_MOVE]= Global.generationCount;
             return _c;
         }
//        if (LxCx.size()<2){
//            //un successfull for this position
//            return c;
//         }
        Lx=LxCx.get(0);
        Cx=LxCx.get(1);
        if(_next<_pos){

            aa[_pos]=new AminoAcid(Cx.x,Cx.y,Cx.z,_c.aminoAcid[_pos].operation);
            for( int count=_pos+1; count<size;count++){
                aa[count]=new AminoAcid(_c.aminoAcid[count-1]);
                //Global.occupiedNode.add(new Node3d(aa[count].x, aa[count].y, aa[count].z));
            }
            aa[_next]=new AminoAcid(Lx.x,Lx.y,Lx.z,_c.aminoAcid[_next].operation);
            for( int count=_next-1; count>=0;count--){
                aa[count]=new AminoAcid(_c.aminoAcid[count+1]);
                //Global.occupiedNode.add(new Node3d(aa[count].x, aa[count].y, aa[count].z));
            }

        }else{
            aa[_pos]=new AminoAcid(Cx.x,Cx.y,Cx.z,_c.aminoAcid[_pos].operation);

            for( int count=_pos-1; count>=0;count--){
                aa[count]=new AminoAcid(_c.aminoAcid[count+1]);
                //Global.occupiedNode.add(new Node3d(aa[count].x, aa[count].y, aa[count].z));
            }
            aa[_next]=new AminoAcid(Lx.x,Lx.y,Lx.z,_c.aminoAcid[_next].operation);
            for( int count=_next+1; count<size;count++){
                aa[count]=new AminoAcid(_c.aminoAcid[count-1]);
                //Global.occupiedNode.add(new Node3d(aa[count].x, aa[count].y, aa[count].z));
            }

        }
         if (Common.checkValidity(aa)==0){
             Global.solExplored+=1;
             RelativeZone rz=new RelativeZone();
             String str=rz.encode(aa);
            //String str=relativeEncode(aa);
//            int fitness=Common.calculateFitness(aa,0);
            int fitness[]={0,0};Common.calculateFitness(aa,fitness);
            Conformation c=new Conformation(aa, fitness, str);
            c.aminoAcid[_pos].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            c.aminoAcid[_next].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            return c;
         }
         else
         {
            _c.aminoAcid[_pos].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            _c.aminoAcid[_next].operation[OperatorType.TILT_MOVE]= Global.generationCount;
            return _c;
         }
    }
//=================tilt move end================================================================

 public  ArrayList<Node3d> getCommonFreeNode(Node3d node1, Node3d node2,Conformation conf){
            ArrayList<Node3d> commonNode=new ArrayList<Node3d>();
            ArrayList<Node3d> freeCommonNode=new ArrayList<Node3d>();
            Node3d[] neighb1=new Node3d[12];
            Node3d[] neighb2=new Node3d[12];

            neighb1=Common.generateTN(node1);
            neighb2=Common.generateTN(node2);

            for (Node3d n1:neighb1){
                for (Node3d n2:neighb2){
                    if(n1.x!=n2.x || n1.y!=n2.y ||n1.z!=n2.z){
                        continue;
                    }
                    commonNode.add(n2);
                }
            }
            for(Node3d n:commonNode){
                boolean exist=false;
                for (AminoAcid a:conf.aminoAcid){
                    if (n.x==a.x && n.y==a.y && n.z==a.z){
                        exist=true;
                        break;
                    }
                }
                if (!exist){
                    freeCommonNode.add(n);
                }
            }


            return freeCommonNode;
        }
     public  ArrayList<Node3d> getFreeNode(Node3d node1,Conformation conf){
        //ArrayList<Node3d> commonNode=new ArrayList<Node3d>();
        ArrayList<Node3d> freeNode=new ArrayList<Node3d>();
        Node3d[] neighb1=new Node3d[12];

        neighb1=Common.generateTN(node1);

        for(Node3d n:neighb1){
            boolean exist=false;
            for (AminoAcid a:conf.aminoAcid){
                if (n.x==a.x && n.y==a.y && n.z==a.z){
                    exist=true;
                    break;
                }
            }
            if (!exist){
                freeNode.add(n);
            }
        }
        return freeNode;
    }


private boolean isFound(Conformation _c,ArrayList<Conformation> neighbor){
    for (Conformation c:neighbor){
        if(c.encodedString.equalsIgnoreCase(_c.encodedString))
            return true;
    }
    return false;
}
}

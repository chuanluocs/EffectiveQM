#ifndef GATE_H_INCLUDED
#define GATE_H_INCLUDED

//CNOT gate
struct gate{
    int control;
    int target;
    int weight; //order in original circuit,
                //0~inf: original gate number,
                //-1:SWAP gate, 
                //-2~-inf: Bridge gate,-x means this Bridge replace x-2th original gate
    // int depth;
    // int type; //0~inf: original gate number, -1:SWAP gate, -2~-inf: Bridge gate,-x means this Bridge replace x-2th original gate
    gate(int c=-1,int t=-1,int w=-1){
        control = c;
        target = t;
        weight = w;
        // depth = d;
        // type = type;
    }
    gate(){
        control = -1;
        target = -1;
        weight = -1;
    }


    bool operator ==(const gate & otherGate) const{
        return control == otherGate.control && target == otherGate.target  && weight == otherGate.weight ;
    }
    bool operator <(const gate & otherGate) const{
        if(control == otherGate.control){
            if(target == otherGate.target){
                return weight < otherGate.weight;
            }
            return target < otherGate.target;
        }
        return control < otherGate.control;
    }
};


#endif // GATE_H_INCLUDED

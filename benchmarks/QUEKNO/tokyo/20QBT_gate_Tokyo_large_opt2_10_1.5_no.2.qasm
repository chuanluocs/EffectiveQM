OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
h q[9];
h q[11];
cx q[12],q[1];
h q[11];
cx q[5],q[2];
h q[16];
h q[1];
cx q[10],q[2];
h q[19];
h q[2];
h q[18];
cx q[10],q[2];
cx q[10],q[8];
h q[5];
h q[1];
h q[8];
h q[2];
cx q[8],q[1];
h q[6];
cx q[11],q[17];
h q[12];
h q[2];
h q[19];
cx q[16],q[0];
cx q[9],q[18];
cx q[12],q[1];
h q[10];
cx q[6],q[19];
cx q[8],q[11];
h q[1];
h q[18];
cx q[5],q[8];
h q[6];
h q[13];
cx q[2],q[1];
cx q[6],q[15];
h q[9];
h q[17];
h q[1];
h q[1];
cx q[5],q[10];
cx q[4],q[9];
h q[13];
h q[17];
h q[17];
cx q[13],q[17];
cx q[5],q[10];
h q[13];
cx q[0],q[6];
h q[13];
cx q[18],q[6];
h q[10];
h q[9];
h q[9];
h q[17];
cx q[0],q[9];
h q[8];
cx q[9],q[6];
cx q[2],q[10];
cx q[3],q[2];
h q[19];
h q[9];
h q[9];
h q[19];
cx q[5],q[1];
h q[5];
h q[9];
cx q[10],q[11];
h q[6];
cx q[2],q[11];
h q[17];
cx q[0],q[6];
cx q[19],q[13];
h q[13];
h q[12];
h q[5];
h q[0];
cx q[2],q[10];
h q[11];
h q[17];
cx q[9],q[6];
cx q[17],q[7];
h q[5];
h q[19];
cx q[4],q[12];
h q[19];
cx q[12],q[1];
h q[7];
h q[19];
cx q[12],q[9];
h q[18];
h q[12];
cx q[4],q[5];
cx q[9],q[5];
h q[16];
cx q[12],q[5];
h q[6];
h q[18];
cx q[0],q[18];
h q[19];
h q[16];
cx q[15],q[19];
h q[15];
h q[5];
cx q[0],q[6];
cx q[15],q[1];
h q[15];
h q[12];
cx q[16],q[0];
h q[18];
cx q[9],q[5];
h q[4];
h q[4];
cx q[0],q[18];
h q[3];
cx q[0],q[5];
h q[3];
h q[0];
cx q[9],q[6];
h q[8];
h q[16];
h q[19];
h q[2];
cx q[12],q[2];
h q[19];
cx q[4],q[12];
cx q[6],q[15];
h q[4];
h q[12];
cx q[3],q[15];
h q[15];
h q[13];
h q[8];
h q[9];
h q[0];
cx q[19],q[13];
h q[4];
cx q[7],q[8];
h q[13];
cx q[16],q[4];
cx q[16],q[4];
cx q[6],q[19];
h q[6];
cx q[5],q[19];
h q[15];
h q[12];
cx q[10],q[2];
cx q[15],q[13];
h q[5];
h q[5];
h q[16];
cx q[7],q[8];
cx q[1],q[7];
h q[16];
h q[13];
h q[5];
cx q[16],q[9];
cx q[7],q[8];
cx q[3],q[1];
h q[2];
h q[10];
h q[2];
cx q[10],q[3];
h q[12];
cx q[12],q[5];
h q[12];
h q[18];
cx q[12],q[5];
h q[19];
h q[12];
h q[19];
h q[9];
h q[15];
h q[18];
cx q[18],q[6];
cx q[18],q[6];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
h q[11];
h q[13];
h q[15];
h q[10];
cx q[10],q[2];
cx q[11],q[13];
cx q[15],q[6];
h q[8];
cx q[2],q[0];
h q[15];
h q[14];
cx q[0],q[4];
cx q[10],q[2];
cx q[12],q[8];
h q[0];
h q[6];
h q[4];
h q[13];
h q[11];
h q[6];
cx q[10],q[2];
h q[2];
cx q[11],q[14];
h q[6];
h q[2];
h q[10];
cx q[15],q[0];
cx q[6],q[2];
h q[15];
cx q[15],q[0];
cx q[15],q[10];
h q[15];
h q[3];
h q[3];
h q[6];
cx q[14],q[3];
h q[6];
cx q[6],q[0];
h q[10];
h q[2];
h q[6];
h q[10];
cx q[15],q[4];
h q[0];
cx q[15],q[2];
cx q[10],q[6];
h q[15];
cx q[15],q[2];
h q[4];
h q[6];
cx q[11],q[12];
cx q[0],q[2];
h q[12];
h q[15];
h q[11];
cx q[6],q[15];
h q[15];
h q[2];
cx q[6],q[15];
h q[6];
cx q[15],q[9];
h q[15];
h q[6];
h q[5];
cx q[15],q[0];
cx q[6],q[2];
h q[0];
h q[8];
cx q[2],q[9];
h q[12];
h q[9];
h q[0];
h q[11];
cx q[14],q[3];
h q[11];
cx q[6],q[2];
h q[5];
cx q[11],q[12];
h q[14];
cx q[8],q[5];
h q[13];
cx q[0],q[2];
h q[0];
h q[2];
cx q[0],q[13];
h q[0];
h q[2];
cx q[0],q[2];
h q[14];
h q[2];
cx q[2],q[4];
cx q[13],q[2];
cx q[12],q[14];
h q[4];
h q[14];
h q[13];
h q[9];
cx q[9],q[4];
h q[4];
cx q[2],q[8];
h q[13];
cx q[9],q[4];
h q[12];
h q[6];
cx q[8],q[9];
cx q[13],q[9];
h q[15];
h q[6];
h q[13];
h q[8];
h q[15];
h q[5];
cx q[10],q[0];
cx q[0],q[15];
h q[1];
h q[8];
cx q[5],q[1];
h q[13];
h q[9];
cx q[15],q[6];
cx q[5],q[1];
h q[7];
h q[12];
cx q[13],q[1];
h q[7];
h q[12];
cx q[9],q[5];
h q[13];
cx q[5],q[7];
cx q[8],q[1];
h q[10];
h q[1];
h q[10];
cx q[8],q[9];
h q[1];
cx q[13],q[1];
h q[5];
h q[14];
h q[8];
cx q[13],q[1];
cx q[10],q[8];
h q[7];
h q[10];
cx q[12],q[14];
cx q[10],q[13];
h q[3];
h q[11];
cx q[13],q[9];
h q[0];
h q[13];
cx q[11],q[14];
h q[3];
cx q[9],q[3];
cx q[9],q[3];
h q[3];
h q[8];
h q[8];
cx q[0],q[11];
cx q[0],q[11];
h q[9];
h q[3];
h q[9];
cx q[8],q[10];
h q[13];
h q[9];
cx q[3],q[7];
h q[10];
h q[13];
cx q[10],q[9];
h q[4];
cx q[15],q[10];
cx q[10],q[9];
h q[4];
cx q[2],q[4];
h q[15];
h q[9];
h q[9];
cx q[10],q[13];
cx q[10],q[13];
h q[3];
h q[4];
h q[4];
h q[3];
h q[4];
cx q[1],q[4];
cx q[15],q[1];
h q[2];
h q[1];
h q[4];
h q[9];
h q[6];
cx q[13],q[4];
h q[1];
cx q[2],q[9];
cx q[6],q[7];
cx q[12],q[13];
h q[7];
h q[12];
cx q[2],q[9];
cx q[12],q[14];
cx q[1],q[2];
cx q[11],q[14];
h q[4];
h q[4];
cx q[13],q[4];
h q[7];
h q[6];
h q[14];
h q[7];
h q[12];
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[13],q[1];
cx q[2],q[13];
h q[13];
h q[14];
h q[1];
cx q[10],q[2];
h q[14];
cx q[5],q[7];
cx q[7],q[0];
h q[7];
h q[15];
cx q[13],q[1];
h q[5];
h q[13];
h q[1];
cx q[15],q[5];
h q[15];
h q[15];
cx q[5],q[14];
h q[7];
cx q[9],q[3];
h q[13];
h q[9];
h q[9];
cx q[5],q[11];
cx q[0],q[6];
h q[1];
h q[3];
h q[12];
cx q[14],q[10];
h q[11];
cx q[0],q[6];
cx q[2],q[13];
h q[2];
cx q[15],q[5];
h q[5];
h q[1];
h q[2];
h q[15];
cx q[6],q[9];
cx q[2],q[12];
h q[5];
h q[9];
h q[1];
cx q[13],q[1];

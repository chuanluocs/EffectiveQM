OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
h q[7];
cx q[15],q[14];
h q[11];
cx q[4],q[6];
h q[8];
h q[4];
cx q[2],q[13];
cx q[4],q[6];
h q[1];
cx q[15],q[8];
cx q[1],q[5];
h q[11];
h q[9];
h q[4];
cx q[16],q[7];
h q[11];
cx q[17],q[19];
cx q[9],q[1];
h q[2];
h q[8];
h q[4];
h q[6];
h q[2];
h q[8];
h q[2];
cx q[7],q[1];
cx q[11],q[8];
h q[16];
h q[7];
h q[1];
cx q[8],q[9];
h q[2];
cx q[17],q[19];

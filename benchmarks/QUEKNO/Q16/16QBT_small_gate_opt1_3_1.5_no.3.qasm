OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[13],q[0];
cx q[5],q[15];
h q[5];
cx q[11],q[8];
h q[13];
h q[13];
h q[0];
h q[15];
h q[11];
cx q[13],q[0];
h q[0];
cx q[13],q[6];
h q[13];
h q[4];
h q[9];
cx q[14],q[8];
h q[5];
cx q[15],q[4];
cx q[5],q[4];
h q[15];
h q[9];
h q[9];
h q[4];
h q[15];
cx q[7],q[9];
cx q[15],q[4];
h q[15];
h q[15];
h q[14];
cx q[8],q[14];
cx q[5],q[15];
h q[13];
h q[8];
h q[12];
cx q[8],q[14];
cx q[12],q[13];
h q[4];
cx q[1],q[2];
cx q[13],q[0];
h q[4];
h q[2];
h q[1];
h q[0];
h q[6];
h q[5];
cx q[13],q[2];
h q[2];
cx q[6],q[0];
cx q[5],q[4];
cx q[1],q[2];
h q[6];

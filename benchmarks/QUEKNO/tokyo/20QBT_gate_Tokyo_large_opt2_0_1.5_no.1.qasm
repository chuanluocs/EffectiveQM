OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[6],q[15];
h q[0];
h q[19];
h q[19];
h q[6];
cx q[13],q[4];
cx q[7],q[0];
h q[12];
cx q[14],q[6];
h q[14];
h q[3];
h q[4];
h q[13];
h q[14];
h q[11];
cx q[2],q[15];
h q[12];
cx q[16],q[18];
cx q[17],q[11];
cx q[12],q[13];
h q[0];
cx q[7],q[19];
h q[2];
h q[19];
h q[15];
h q[13];
cx q[12],q[13];
h q[0];
cx q[13],q[4];
cx q[3],q[19];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
h q[7];
cx q[7],q[12];
h q[15];
h q[15];
cx q[14],q[15];
h q[12];
cx q[14],q[15];
h q[7];
h q[15];
h q[7];
cx q[10],q[15];
h q[7];
cx q[14],q[10];
h q[15];
h q[10];
cx q[3],q[6];
h q[0];
cx q[6],q[7];
h q[3];
h q[3];
cx q[10],q[15];
cx q[5],q[0];
h q[10];
h q[15];
cx q[5],q[11];
h q[4];
h q[3];
h q[4];
h q[8];
h q[15];
h q[7];
h q[10];
cx q[10],q[4];
cx q[8],q[5];
cx q[15],q[3];
h q[3];
h q[8];
cx q[11],q[0];
cx q[8],q[5];
h q[12];
cx q[7],q[12];

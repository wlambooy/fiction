// Benchmark "ex" written by ABC on Fri Jul  4 15:14:14 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8;
  assign new_n6 = a & ~b;
  assign new_n7 = a & ~new_n6;
  assign new_n8 = ~d & ~new_n7;
  assign F0 = ~c & new_n8;
endmodule



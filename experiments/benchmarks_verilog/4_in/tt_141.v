// Benchmark "ex" written by ABC on Fri Jul  4 15:14:13 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7;
  assign new_n6 = ~a & ~b;
  assign new_n7 = ~c & new_n6;
  assign F0 = ~d & new_n7;
endmodule



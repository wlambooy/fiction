// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b,
    F0  );
  input  a, b;
  output F0;
  wire new_n4, new_n5;
  assign new_n4 = a & ~b;
  assign new_n5 = ~a & ~b;
  assign F0 = new_n4 | new_n5;
endmodule



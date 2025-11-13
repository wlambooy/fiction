// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b, c,
    F0  );
  input  a, b, c;
  output F0;
  wire new_n5;
  assign new_n5 = ~a & ~b;
  assign F0 = ~c & new_n5;
endmodule



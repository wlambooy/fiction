// Benchmark "ex" written by ABC on Fri Jul  4 15:14:12 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6;
  assign new_n6 = ~c & ~d;
  assign F0 = ~b & new_n6;
endmodule



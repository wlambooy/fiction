// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b, c,
    F0  );
  input  a, b, c;
  output F0;
  wire new_n5, new_n6;
  assign new_n5 = ~b & c;
  assign new_n6 = b & ~c;
  assign F0 = new_n5 | new_n6;
endmodule



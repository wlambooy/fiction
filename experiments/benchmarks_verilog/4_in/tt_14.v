// Benchmark "ex" written by ABC on Fri Jul  4 15:14:13 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7;
  assign new_n6 = ~c & d;
  assign new_n7 = c & ~d;
  assign F0 = new_n6 | new_n7;
endmodule



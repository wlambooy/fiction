module top (
    A,     // Input for FO1 and FO2
    B,     // Input for FO1 and FO2
    AND_out,  // Output from AND gate
    OR_out    // Output from OR gate
);

input A, B;
output AND_out, OR_out;

// FO1 Logic
wire FO1_inv;  // Inverted signal from FO1
wire FO1_ho;   // Direct signal from FO1

// FO2 Logic
wire FO2_buf;  // Buffered signal from FO2
wire FO2_ho;   // Direct signal from FO2

// FO1 logic
assign FO1_inv = ~A;   // Invert input A for FO1 (INV)
assign FO1_ho = A;     // Direct connection from FO1

// FO2 logic
assign FO2_buf = B;    // Buffered input B for FO2 (BUF)
assign FO2_ho = B;     // Direct connection from FO2

// AND gate logic
assign AND_out = FO1_buf & FO1_ho;  // FO1 → BUF → AND, FO1 → HO → AND

// OR gate logic
assign OR_out = FO1_inv | FO2_ho;   // FO1 → INV → OR, FO2 → HO → OR

endmodule
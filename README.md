# HAlign3.5

HAlign3.5 is a high-performance multiple sequence alignment software based on the star alignment strategy, designed for efficiently aligning large numbers of sequences. Compared to its predecessor HAlign3, HAlign3.5 further enhances the ability to handle long sequences and large-scale datasets, enabling fast and efficient alignment on standard computing devices.

## Background
[HAlign3](https://github.com/malabz/HAlign-3) was implemented in Java and was capable of efficiently aligning ultra-large sets of similar DNA/RNA sequences, but had limitations when dealing with long sequences and very large datasets. To address these issues, HAlign3.5 was reimplemented in C++ and incorporates the Burrows-Wheeler Transform (BWT) and wavefront alignment algorithm.

### Key Improvements
- **Algorithm Optimization**: Replaced the original suffix tree with BWT for more efficient indexing and searching.
- **Memory and Speed Optimization**: Introduced the wavefront alignment algorithm to reduce memory usage and improve alignment speed, especially for long sequences.


## Compilation
HAlign3.5 is written in C++ and can be compiled using the `make` tool.

```bash
make
```

After compilation, an executable file named `halign3.5` will be generated.

## Usage
```bash
./halign3.5 Input_file Output_file [-r/--reference val] [-t/--threads val] [-sa/--sa val] [-h/--help]
```

### Parameter Description
- `Input_file`: Path to the input file or folder (please use `.fasta` as the file suffix or input folder).
- `Output_file`: Path to the output file (please use `.fasta` as the file suffix).
- `-r/--reference`: Reference sequence name (please remove all whitespace), default is the longest sequence.
- `-t/--threads`: Number of threads to use, default is 1.
- `-sa/--sa`: Global `sa` threshold, default is 15.
- `-h/--help`: Show help information.

## Example
Here is a simple example of using HAlign3.5 for multiple sequence alignment:

```bash
./halign3.5 input.fasta output.fasta -t 4
```

This command will use 4 threads to align the sequences in the `input.fasta` file, and the result will be saved in the `output.fasta` file.

## Reference
- [HAlign Official Website](http://lab.malab.cn/soft/halign/)

## License
HAlign3.5 is developed by the Malab team under the [MIT License](https://github.com/metaphysicser/HAlign3.5/blob/main/LICENSE).
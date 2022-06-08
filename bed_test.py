import unittest
import io
import bed

example_bed = """
test_chrom1 20 30 name1 0 +
test_chrom2 30 40 name2 0 -
test_chrom3 40 50 name3 0 .
""".lstrip()  # remove new line at start


class TestBED(unittest.TestCase):
    def test_bed_reader(self):
        example_file = io.StringIO(example_bed)
        bed_reader = bed.Reader(example_file)

        for genomic_region in bed_reader:
            print(genomic_region)

    def test_bed_writer(self):
        example_file = io.StringIO("")

        bed_writer = bed.Writer(example_file)
        bed_writer.write_genomic_region(
            bed.GenomicRegion(chrom="test", start=5, end=10)
        )

        example_file.seek(0)
        print(example_file.read())


if __name__ == "__main__":
    unittest.main()

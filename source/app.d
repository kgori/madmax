import bio.bam.pileup;
import bio.bam.reader;
import std.algorithm: copy, count, filter, isSorted, map, sort, topN, uniq;
import std.array;
import std.conv: to;
import std.file: exists;
import std.getopt;
import std.math: fabs;
import std.range;
import std.stdio;

static string usage = "MADmax - scan for amplified read counts using MAD\nUsage: madmax -b|--bamfile [PATH] -c|--mad-const [FLOAT] -w|--window-size [INT] -r|min-run [INT]";

struct Position {
    ulong depth;
    ulong pos;
    ulong ref_id;
}

struct PositionRange {
    ulong ref_id;
    ulong start;
    ulong end;
    ulong depthsum;
    ulong n;
}

bool read_filter(ReadType)(ReadType read) {
    return !(read.is_duplicate() || read.failed_quality_control() || read.current_base == '-');
}

auto deviation_fn(ulong median) {
    return (ulong val){ return median > val ? median - val : val - median; };
}

/++
Reads a pileup region (from a Range of columns) into a buffer (which
        is an array of Positions)
+/
void read_region(PileupColumns)(PileupColumns columns, ref Position[] buffer, int buflen) {
    foreach (pos; buffer) pos = Position();
    foreach (i, column; enumerate(columns.take(buflen))) {
        if (i >= buffer.length) break;
        ulong depth = column.reads.count!(read => read_filter(read));
        ulong pos = column.position;
        ulong ref_id = column.ref_id;
        buffer[i] = Position(depth, pos, ref_id);
    }
}


double find_median(Numeric)(ref Numeric[] array) {
    assert(isSorted(array));
    if (array.length % 2 == 0) {
        return (to!double(array[(array.length/2) - 1]) + to!double(array[array.length/2])) / 2;
    }
    else {
        return to!double(array[(array.length/2)]);
    }
}

unittest {
    ulong[] even = [1,4,5,17,18,23,29,50,61,62];
    assert(find_median(even) == (18. + 23.) / 2);
    ulong[] odd = [1,4,5,17,18,23,29,50,61,62,75];
    assert(find_median(odd) == 23.);
}


Array find_outliers(Array)(ref Array buf1, ref Array buf2, double constval) {
    ulong start = buf1[0].pos;
    ulong end = buf2[$-1].pos;
    auto depths = chain(buf1, buf2)
        .filter!"a.pos > 0"
        .map!(item => item.depth)
        .array
        .sort
        .array;
    double median_depth = find_median(depths);

    auto deviations = chain(buf1, buf2)
        .filter!(a => (a.pos > 0) && (a.depth > median_depth))
        .map!(item => item.depth > median_depth ? item.depth - median_depth : median_depth - item.depth)
        .array
        .sort
        .array;
    double mad = find_median(deviations);

    // writefln("DEBUG:region %d-%d: N = %d; median depth = %f; mad = %f", start, end, deviations.length, median_depth, mad);

    return chain(buf1, buf2)
        .filter!(item => (item.depth - median_depth) > mad*constval)
        .array;
}

/++
Returns true if b is consecutive to a
+/
bool consecutive(Position a, Position b) {
    return (a.ref_id == b.ref_id) && (a.pos + 1 == b.pos);
}

void write_outliers(Array)(Array outliers1, Array outliers2, int minRunLength, ref BamReader bam) {
    auto prev = Position();
    auto range = PositionRange(prev.ref_id, prev.pos, prev.pos, 0, 1);
    foreach (pos; chain(outliers1, outliers2).sort!"a.pos < b.pos".uniq) {
        if (consecutive(prev, pos)) {
            // extend current range
            range.end++;
            range.n++;
            range.depthsum += pos.depth;
        }
        else {
            // emit current range
            if ((range.end - range.start + 1) > minRunLength && range.depthsum > 0) {
                writefln("%s:%d-%d\t%f",
                         bam.reference_sequences[range.ref_id].name, range.start+1, range.end+1,
                         to!double(range.depthsum) / range.n);
            }
            range.ref_id = pos.ref_id;
            range.start = pos.pos;
            range.end = pos.pos;
        }
        prev = pos;
    }
    // emit any remaining range
    if ((range.end - range.start + 1) > minRunLength && range.depthsum > 0) {
        writefln("%s:%d-%d\t%f",
                 bam.reference_sequences[range.ref_id].name, range.start+1, range.end+1,
                 to!double(range.depthsum) / range.n);
    }
}

void main(string[] argv)
{
    string bamfile;
    int windowSize = 10000;
    int minRun = 5;
    double madConst = 1.4826;

    try {
        auto args = getopt(
                argv,
                std.getopt.config.required, "bamfile|b", "Path to bamfile", &bamfile,
                "window-size|w", "Size of window used to compute MAD [Default 10000]", &windowSize,
                "min-run|r", "Minimum number of consecutive position above MAD required to call region [Default: 5]", &minRun,
                "mad-const|c", "Outliers are defined as those with MAD above (mad-const * MAD) [Default: 1.4826]", &madConst,
                );

        if (args.helpWanted) {
            defaultGetoptPrinter(usage, args.options);
            return;
        }
    }
    catch (GetOptException) {
        writeln(usage);
    }

    int BUFLEN = windowSize/2;

    if (!exists(bamfile)) {
        writefln("File %s does not exist. Exiting", bamfile);
        return;
    }

    auto bam = new BamReader(bamfile);

    if (!bam.has_index()) {
        bam.createIndex();
    }

    // TODO: Use running median using an Indexable Skiplist: https://code.activestate.com/recipes/576930/
    // for now, this is the dumb version
    foreach (ref reference; bam.reference_sequences) {
        auto buffer1 = new Position[BUFLEN];
        auto buffer2 = new Position[BUFLEN];
        auto buffer3 = new Position[BUFLEN];
        auto pileup = makePileup(bam[reference.name][]);
        bool initialised = false;

        if (!pileup.empty && !initialised) {
            read_region(refRange(&pileup), buffer1, BUFLEN);
            read_region(refRange(&pileup), buffer2, BUFLEN);
            read_region(refRange(&pileup), buffer3, BUFLEN);
            auto outliers1 = find_outliers(buffer1, buffer2, madConst);
            auto outliers2 = find_outliers(buffer2, buffer3, madConst);
            write_outliers(outliers1, outliers2, minRun, bam);
            initialised = true;
        }

        while (!pileup.empty) {
            buffer1 = buffer2;
            buffer2 = buffer3;
            read_region(refRange(&pileup), buffer3, BUFLEN);
            auto outliers1 = find_outliers(buffer1, buffer2, madConst);
            auto outliers2 = find_outliers(buffer2, buffer3, madConst);
            write_outliers(outliers1, outliers2, minRun, bam);
        }
    }
}

import idautils
import ida_ida
import ida_search
import ida_bytes
import ida_kernwin

import typing
import functools
import itertools


class ScoreUnit(object):
    W1 = 0
    W2 = 1
    W3 = 2
    W4 = 3
    W1_MASK = 1 << W1
    W2_MASK = 1 << W2
    W3_MASK = 1 << W3
    W4_MASK = (1 << 8) - (1 << W4)

    # def __init__(self, score: int, x: int, y: int, p, target_gap: int = 0, match: int = 0, diff_gap: int = 0,
    #              diff_combine: int = 0):

    # we only choose one path of same score, so every unit has one parent
    def __init__(self, score: int, x: int, y: int, p):
        self.score = score
        # self.W = (target_gap << self.W1) + (match << self.W2) + (diff_gap << self.W3) + (diff_combine << self.W4)
        self.x = x
        self.y = y
        self.p = p
        if p is not None:
            p.child.append(self)
            # p.add_child(self)
        self.child = []

    # def position(self) -> tuple:
    #     return self.x, self.y
    #
    # def get_parent(self):
    #     return self.p
    #
    # def get_child(self):
    #     return self.child
    #
    # def add_child(self, child):
    #     self.child.append(child)


class MatrixException(Exception): pass


class Matrix(object):
    class MatrixN(object):
        def __init__(self, n):
            self.n = n
            self.y: typing.Dict[int, ScoreUnit] = {}

        def __getitem__(self, item: int):
            return self.y[item]

        def __setitem__(self, key: int, value: ScoreUnit):
            if key >= self.n:
                raise MatrixException("out of scope")
            self.y[key] = value

        def __delitem__(self, key: int):
            del self.y[key]

    def __init__(self, m: int, n: int):
        self.m = m
        self.n = n
        self.m_array: typing.List = [self.MatrixN(self.n) for i in range(m)]

    def __getitem__(self, item):
        return self.m_array[item]

    def __setitem__(self, key, value):
        self.m_array[key] = value

    def __delitem__(self, key):
        del self.m_array[key]


def clean_leaf(root: ScoreUnit, now_y: int):
    stack_chunks = []

    child_size = len(root.child)
    current_i = 0
    current_unit: ScoreUnit = root

    while True:
        if current_i < child_size:
            p = current_unit.child[current_i]
            stack_chunks.append((current_unit, current_i, child_size))
            child_size = len(p.child)
            current_i = 0
            current_unit = p
        else:
            if len(stack_chunks) == 0:
                break  # travel finish
            # if we here, this mean no child or child iter end, pop a chunk
            p = current_unit
            current_unit, current_i, child_size = stack_chunks.pop(-1)
            if len(p.child) == 0 and p.y != now_y and p.x != 0 and p.y != 0:
                current_unit.child.pop(current_i)
                child_size -= 1
            else:
                current_i += 1


# def clean_leaf(root: ScoreUnit, now_y: int):
#     for node in root.child:
#         clean_leaf(node, now_y)
#
#     root.child = filter(lambda n: len(n.child) == 0 and n.y != now_y, root.child)


def get_align_list(matrix: Matrix, target_seq: list, diff_seq: list):
    m_target_seq = []
    m_diff_seq = []
    n: ScoreUnit = matrix[matrix.m - 1][matrix.n - 1]
    while n.x != 0 and n.y != 0:
        p: ScoreUnit = n.p
        x_list = list(zip(target_seq[p.x: n.x], range(p.x, n.x)))
        x_list.reverse()
        y_list = list(zip(diff_seq[p.y: n.y], range(p.y, n.y)))
        y_list.reverse()
        for x_v, y_v in itertools.zip_longest(x_list, y_list, fillvalue=(0, 0)):
            m_target_seq.append(x_v)
            m_diff_seq.append(y_v)
        n = p

    m_target_seq.reverse()
    m_diff_seq.reverse()
    return m_target_seq, m_diff_seq


def seq_align(target_seq, diff_seq):
    x_gap = -3
    y_gap = -7
    s_match = 10
    s_mismatch = -8
    s_combine = -5

    m = len(target_seq)
    n = len(diff_seq)
    matrix = Matrix(m, n)

    # init matrix
    matrix[0][0] = ScoreUnit(0, 0, 0, None)

    prev: ScoreUnit = matrix[0][0]
    for i in range(1, m):
        matrix[i][0] = prev = ScoreUnit(x_gap * i, i, 0, prev)

    prev = matrix[0][0]
    for j in range(1, n):
        matrix[0][j] = prev = ScoreUnit(y_gap * j, 0, j, prev)

    for y in range(1, n):
        print(y)
        for x in range(1, m):
            s1: ScoreUnit = matrix[x - 1][y]
            score1 = s1.score + x_gap, s1
            s2: ScoreUnit = matrix[x][y - 1]
            score2 = s2.score + y_gap, s2
            s3: ScoreUnit = matrix[x - 1][y - 1]
            if target_seq[x - 1] == diff_seq[y - 1]:
                score3 = s3.score + s_match, s3
            else:
                diff = diff_seq[y - 1]
                val = target_seq[x - 1]
                i = 1
                while i < 3 and x - 1 - i >= 0:  # TODO: s_match =  -10 * s_combine
                    if val >= diff:
                        break
                    val += target_seq[x - 1 - i]
                    i += 1

                if val == diff:
                    s4: ScoreUnit = matrix[x - 1 - i][y - 1]
                    score3 = s4.score + s_match + i * s_combine, s4
                else:
                    score3 = s3.score + s_mismatch, s3
            score, s = max(score1, score2, score3, key=lambda t: t[0])
            matrix[x][y] = ScoreUnit(score, x, y, s)
        clean_leaf(matrix[0][0], y)
    return get_align_list(matrix, target_seq, diff_seq)


def adjacent_difference(seq: typing.Iterable, op: typing.Callable = None) -> typing.Generator[int, None, None]:
    acc = next(iter(seq))
    for val in seq:
        if op is None:
            yield val - acc
        else:
            yield op(val - acc)
        acc = val


# first get all symbol from idb
def get_symbol_addr() -> typing.Generator[int, None, None]:
    for s in idautils.Strings():
        yield s.ea


def get_all_point_data() -> typing.Generator[int, None, None]:
    ea = ida_ida.inf_get_min_ea()
    while ea < ida_ida.inf_get_max_ea():
        ea = ida_search.find_data(ea, ida_search.SEARCH_NEXT | ida_search.SEARCH_DOWN)
        if ida_bytes.get_item_size(ea) == 4:
            yield ida_bytes.get_wide_dword(ea)


def choose_point_seq(point_diff_seq: list) -> typing.Tuple[typing.List[int], int]:
    # point_diff_seq = [i for i in adjacent_difference(get_all_point_data())]
    seq_size = len(point_diff_seq)
    if seq_size <= 10000:
        min_idx = 0
        point_seq = point_diff_seq
    else:
        # match to large sequence will be slow and use many memory, limit to 10000
        # choose the most "close" data
        min_avg = sum(itertools.islice(point_diff_seq, 10000))  # user iter to save some mem
        min_idx = 0

        prev_avg = min_avg
        for i in range(1, seq_size - 9999):
            now_avg = prev_avg + point_diff_seq[i + 9999] - point_diff_seq[i - 1]
            if now_avg < min_avg:
                min_avg = now_avg
                min_idx = i
            prev_avg = now_avg

        point_seq = point_diff_seq[min_idx: min_idx + 10000]
    return point_seq, min_idx


def main():
    sym_list = list(get_symbol_addr())
    target_seq = adjacent_difference(sym_list)

    all_point = list(set(get_all_point_data()))
    all_point.sort()

    diff_all_point = list(adjacent_difference(all_point))
    diff_seq, idx = choose_point_seq(diff_all_point)

    align_target, align_diff = seq_align(list(target_seq), diff_seq[:200])

    for x, y in align_target[:1000]:
        print("%03d" % x, end=', ')
    print("\n")
    for x, y in align_diff[:1000]:
        print("%03d" % x, end=', ')
    print("\n")

    base_offs = {}
    for point_diff, sym_diff in filter(lambda x: x[0][0] != 0 and x[1][0] != 0, zip(align_diff, align_target)):
        print(point_diff, sym_diff)
        point_idx = point_diff[1] - 1
        sym_idx = sym_diff[1] - 1
        sym_addr = sym_list[sym_idx + 1]
        point_addr = all_point[idx + point_idx + 1]
        print(point_addr, all_point[idx + point_idx], sym_addr, sym_list[sym_idx])
        base_off = point_addr - sym_addr
        if base_off in base_offs:
            base_offs[base_off] += 1
        else:
            base_offs[base_off] = 1
    max_count = 0
    max_count_off = 0
    print(base_offs)
    for off, count in base_offs.items():
        if count > max_count:
            max_count = count
            max_count_off = off

    print(max_count)
    print(max_count_off)


if __name__ == '__main__':
    main()

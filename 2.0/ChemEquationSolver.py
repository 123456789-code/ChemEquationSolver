def add_dict(d1, d2):
    for k, v in d2.items():
        if k in d1:
            d1[k] += v
        else:
            d1[k] = v
    return d1


def times_dict(d, n):
    for i in d.keys():
        d[i] *= n
    return d


class RationalNumber:
    def __init__(self, p=0, q=1):
        self.p, self.q = p, q

    def gcd(self, a, b):
        (a, b) = (a, b) if a >= b else (b, a)
        if b == 0:
            return 1
        else:
            while a % b != 0:
                a, b = b, a % b
            return b

    def simplify(self):
        if self.p == 0:
            self.q = 1
        else:
            gcd = self.gcd(abs(self.p), abs(self.q))
            self.p //= gcd
            self.q //= gcd
            if self.q < 0:
                self.p, self.q = -self.p, -self.q
        return self

    def __add__(self, other):
        return RationalNumber(self.p * other.q + self.q * other.p, self.q * other.q).simplify()

    def __sub__(self, other):
        return RationalNumber(self.p * other.q - self.q * other.p, self.q * other.q).simplify()

    def __mul__(self, other):
        return RationalNumber(self.p * other.p, self.q * other.q).simplify()

    def __truediv__(self, other):
        if other.p == 0:
            raise ZeroDivisionError
        else:
            return RationalNumber(self.p * other.q, self.q * other.p).simplify()

    def __iadd__(self, other):
        return self + other

    def __isub__(self, other):
        return self - other

    def __imul__(self, other):
        return self * other

    def __itruediv__(self, other):
        return self / other

    def __lt__(self, other):
        return (self - other).p < 0

    def __le__(self, other):
        return (self - other).p <= 0

    def __gt__(self, other):
        return (self - other).p > 0

    def __ge__(self, other):
        return (self - other).p >= 0

    def __eq__(self, other):
        return self.p == other.p and self.q == other.q

    def __ne__(self, other):
        return not self == other

    def __float__(self):
        return self.p / self.q

    def __str__(self):
        return str(self.p) + "/" + str(self.q) if self.p != 0 and self.q != 1 else str(self.p)

    def __int__(self):
        if self.p == 0 or self.q == 1:
            return self.p
        else:
            raise RuntimeError("rational_number in int() is not an integer")

    def load(self, p, q=1):
        if q == 0:
            raise ZeroDivisionError
        else:
            self.p, self.q = p, q
            self.simplify()

    def load_string(self, string):
        if "/" in string:
            [p, q] = string.split("/")
            self.load(int(p), int(q))
        else:
            self.load(int(string))


class Matrix:
    def __init__(self, lst):
        self.list = []
        for i in lst:
            lst_cache = []
            for j in i:
                lst_cache.append(RationalNumber(j))
            self.list.append(lst_cache)

    def __str__(self):
        string = ""
        for i in self.list:
            for j in i:
                string += str(j) + ", "
            string += "\n"
        return string

    def del_zero(self):
        lst = []
        for i in range(len(self.list)):
            if self.list[i] == len(self.list[i]) * [RationalNumber(0)]:
                lst.append(i)
        for i in range(len(lst)):
            lst[i] -= i
        for i in lst:
            del self.list[i]
        return self

    def row_echelon_form(self):
        # 行阶梯矩阵
        for i in range(len(self.list[0])):
            if i >= len(self.list):
                break
            if self.list[i][i] == RationalNumber(0):
                for j in range(i + 1, len(self.list)):
                    if self.list[j][i] != RationalNumber(0):
                        self.list[i], self.list[j] = self.list[j], self.list[i]
                        break
                if self.list[i][i] == RationalNumber(0):
                    continue
            for j in range(i + 1, len(self.list)):
                factor = self.list[j][i] / self.list[i][i]
                for k in range(i, len(self.list[0])):
                    self.list[j][k] -= factor * self.list[i][k]
        self.del_zero()

    def solve(self):
        self.row_echelon_form()
        if len(self.list[0]) > len(self.list) + 1:
            return ["Infinity_solution"]
        elif len(self.list[0]) < len(self.list) + 1:
            return ["No_solution"]
        else:
            coefficient = [RationalNumber(1)]
            for i in range(len(self.list)):
                cache = RationalNumber(0)
                for j in range(i + 1):
                    cache -= self.list[-1 - i][-1 - j] * coefficient[j]
                if self.list[-1 - i][-2 - j] == RationalNumber(0):
                    if cache == 0:
                        return ["Infinity_solution"]
                    else:
                        coefficient = [RationalNumber(0)] * (i + 1) + [RationalNumber(1)]
                else:
                    coefficient.append(cache / self.list[-1 - i][-2 - j])
            coefficient.reverse()
            lcm = 1
            for i in coefficient:
                lcm = lcm * i.q // i.gcd(lcm, i.q)
            lcm_ = RationalNumber(-lcm)
            for i in range(len(coefficient)):
                coefficient[i] *= lcm_
            return coefficient


class Thing:
    def __init__(self, string):
        self.string = string
        string = string.replace(" ", "").replace("\n", "")
        if string == "e[-]":
            self.elements = {"e": 0}
            self.electron = -1
        elif "]" in string:
            lst = string[:-1].split("[")
            if lst[1] == "":
                self.electron = 0
            elif lst[1] == "+":
                self.electron = 1
            elif lst[1] == "-":
                self.electron = -1
            else:
                self.electron = int(lst[1][-1] + lst[1][:-1])
            self.elements = self.separate(lst[0].replace(" ", "").replace("\n", ""), 0)[0]
        else:
            self.electron = 0
            self.elements = self.separate(string.replace(" ", "").replace("\n", ""), 0)[0]

    def separate(self, string, position):
        element_dict = {}
        name_cache = ""
        num_cache = ""
        dict_cache = {}
        while position < len(string):
            char = string[position]
            if char == "(":
                if name_cache != "":
                    element_dict = add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
                    name_cache = ""
                    num_cache = ""
                dict_cache, position = self.separate(string, position + 1)
            elif char == ")":
                break
            elif char.isdigit():
                num_cache += char
            elif char.islower():
                name_cache += char
            elif char.isupper:
                if dict_cache == {}:
                    element_dict = add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
                    name_cache = char
                    num_cache = ""
                else:
                    element_dict = add_dict(element_dict,
                                            times_dict(dict_cache, (1 if num_cache == "" else int(num_cache))))
                    dict_cache = {}
                    name_cache = char
                    num_cache = ""
            else:
                raise RuntimeError("非法字符")
            position += 1
        element_dict = add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
        element_dict = add_dict(element_dict, times_dict(dict_cache, (1 if num_cache == "" else int(num_cache))))
        if "" in element_dict:
            del element_dict[""]
        return element_dict, position


class Equation:
    def __init__(self, string):
        self.things = []
        lst = string.split("--")
        for i in lst[0].split(" + ") + lst[1].split(" + "):
            i = i.replace(" ", "").replace("\n", "")
            self.things.append(Thing(i))
        self.element_list = []
        for i in self.things:
            self.element_list += i.elements.keys()
        self.element_list = list(set(self.element_list))
        self.coefficient = []

    def solve(self):
        matrix_cache = []
        for i in self.element_list:
            lst = []
            for j in self.things:
                lst.append(j.elements[i] if i in j.elements else 0)
            matrix_cache.append(lst)
        lst = []
        for j in self.things:
            lst.append(j.electron)
        matrix_cache.append(lst)
        matrix1 = Matrix(matrix_cache).del_zero()
        self.coefficient = matrix1.solve()

    def __str__(self):
        if self.coefficient == ["No_solution"]:
            return "无解"
        elif self.coefficient == ["Infinity_solution"]:
            return "无穷解"
        else:
            left_string, right_string = "", ""
            for i in range(len(self.coefficient)):
                if self.coefficient[i] > RationalNumber(0):
                    if self.coefficient[i] != RationalNumber(1):
                        left_string += str(self.coefficient[i])
                    left_string += self.things[i].string + " + "
                elif self.coefficient[i] < RationalNumber(0):
                    if self.coefficient[i] != RationalNumber(-1):
                        right_string += str(self.coefficient[i])[1:]
                    right_string += self.things[i].string + " + "
            return left_string[:-3] + " == " + right_string[:-3]


a = input()
equation_in = Equation(a)
equation_in.solve()
print(equation_in)

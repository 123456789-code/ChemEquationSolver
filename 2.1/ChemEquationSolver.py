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

    def __abs__(self):
        return RationalNumber(abs(self.p), self.q)

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

    def isinteger(self):
        return self.q == 1

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
    def __init__(self, lst, outcome):
        self.list = lst
        self.outcome = outcome

    def __str__(self):
        string = ""
        for i in range(len(self.list)):
            for j in self.list[i]:
                string += str(j) + ", "
            string += "|| " + str(self.outcome[i]) + "\n"
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
            del self.outcome[i]
        return self

    def row_echelon_form(self):  # 行阶梯矩阵
        for i in range(len(self.list[0])):
            if i >= len(self.list):
                break
            if self.list[i][i] == RationalNumber(0):
                for j in range(i + 1, len(self.list)):
                    if self.list[j][i] != RationalNumber(0):
                        self.list[i], self.list[j] = self.list[j], self.list[i]
                        self.outcome[i], self.outcome[j] = self.outcome[j], self.outcome[i]
                        break
                if self.list[i][i] == RationalNumber(0):
                    continue
            for j in range(i + 1, len(self.list)):
                factor = self.list[j][i] / self.list[i][i]
                for k in range(i, len(self.list[0])):
                    self.list[j][k] -= factor * self.list[i][k]
                self.outcome[j] -= factor * self.outcome[i]
        self.del_zero()

    def solve(self):
        self.row_echelon_form()
        # print(self)
        if len(self.list[0]) > len(self.list):
            return ["Infinity_solution"]
        elif len(self.list[0]) < len(self.list):
            return ["No_solution"]
        else:
            coefficient = [self.outcome[-1] / self.list[-1][-1]]
            # print("coefficient.append:", coefficient[0])
            del self.list[-1]
            del self.outcome[-1]
            for i in range(len(self.list)):
                cache = self.outcome[-1 - i]
                for j in range(i + 1):
                    # print("cache -=", self.list[-1 - i][-1 - j] * coefficient[j])
                    cache -= self.list[-1 - i][-1 - j] * coefficient[j]
                if self.list[-1 - i][-2 - j] == RationalNumber(0):
                    if cache == 0:
                        return ["Infinity_solution"]
                    else:
                        coefficient = [RationalNumber(0)] * (i + 1) + [RationalNumber(1)]
                else:
                    # print("coefficient.append:", cache / self.list[-1 - i][-2 - j])
                    coefficient.append(cache / self.list[-1 - i][-2 - j])
            coefficient.reverse()
            # for i in coefficient:
            #     print(i,end=", ")
            # print("\n")
            return coefficient


class Thing:
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
                    element_dict = self.add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
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
                    element_dict = self.add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
                    name_cache = char
                    num_cache = ""
                else:
                    element_dict = self.add_dict(element_dict,
                                                 self.times_dict(dict_cache,
                                                                 (1 if num_cache == "" else int(num_cache))))
                    dict_cache = {}
                    name_cache = char
                    num_cache = ""
            else:
                raise RuntimeError("非法字符")
            position += 1
        element_dict = self.add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
        element_dict = self.add_dict(element_dict,
                                     self.times_dict(dict_cache, (1 if num_cache == "" else int(num_cache))))
        if "" in element_dict:
            del element_dict[""]
        return element_dict, position


class Equation:
    def __init__(self, string):
        self.things = []
        self.coefficient_in = []
        self.coefficient = []
        lst = string.split("--")
        for i in lst[0].split(" + "):
            i = i.replace(" ", "").replace("\n", "")
            if i[0].isdigit() or i[0] == "-":
                cache = 0
                while i[cache].isdigit() or i[cache] in "/-":
                    cache += 1
                if "/" in i[:cache]:
                    lst_cache = i[:cache].split("/")
                    self.coefficient_in.append(RationalNumber(int(lst_cache[0]), int(lst_cache[1])))
                else:
                    self.coefficient_in.append(RationalNumber(int(i[:cache])))
                self.things.append(Thing(i[cache:]))
            else:
                self.coefficient_in.append(RationalNumber(0))
                self.things.append(Thing(i))
        for i in lst[1].split(" + "):
            i = i.replace(" ", "").replace("\n", "")
            if i[0].isdigit() or i[0] == "-":
                cache = 0
                while i[cache].isdigit() or i[cache] in "/-":
                    cache += 1
                if "/" in i[:cache]:
                    lst_cache = i[:cache].split("/")
                    self.coefficient_in.append(RationalNumber(-int(lst_cache[0]), int(lst_cache[1])))
                else:
                    self.coefficient_in.append(RationalNumber(-int(i[:cache])))
                self.things.append(Thing(i[cache:]))
            else:
                self.coefficient_in.append(RationalNumber(0))
                self.things.append(Thing(i))
        self.element_list = []
        for i in self.things:
            self.element_list += i.elements.keys()
        self.element_list = list(set(self.element_list))

    def solve(self):
        matrix_cache = []
        for i in self.element_list:
            lst = []
            for j in self.things:
                lst.append(RationalNumber(j.elements[i]) if i in j.elements else RationalNumber(0))
            matrix_cache.append(lst)
        lst = []
        for i in self.things:
            lst.append(RationalNumber(i.electron))
        matrix_cache.append(lst)
        outcome_cache = [RationalNumber(0)] * len(matrix_cache)
        flag = True
        for i in range(len(self.coefficient_in)):
            if self.coefficient_in[i] != RationalNumber(0):
                matrix_cache.append([RationalNumber(0)] * i + [RationalNumber(1)] +
                                    [RationalNumber(0)] * (len(matrix_cache[0]) - i - 1))
                outcome_cache.append(self.coefficient_in[i])
                flag = False
        if flag:
            matrix_cache.append([RationalNumber(0)] * (len(matrix_cache[0]) - 1) + [RationalNumber(1)])
            outcome_cache.append(RationalNumber(1))
        matrix1 = Matrix(matrix_cache, outcome_cache).del_zero()
        self.coefficient = matrix1.solve()
        if flag:
            lcm = 1
            for i in self.coefficient:
                lcm = lcm * i.q // i.gcd(lcm, i.q)
            lcm_ = RationalNumber(-lcm)
            for i in range(len(self.coefficient)):
                self.coefficient[i] *= lcm_
        if not flag:
            for i in range(len(self.coefficient)):
                if self.coefficient_in[i] != RationalNumber(0) and self.coefficient_in[i] != self.coefficient[i]:
                    self.coefficient = ["No_solution"]
                    break

    def __str__(self):
        # for i in self.coefficient:
        #     print(i, end=", ")
        # print("\n")
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
                        if not self.coefficient[i].isinteger():
                            left_string += " "
                    left_string += self.things[i].string + " + "
                elif self.coefficient[i] < RationalNumber(0):
                    if self.coefficient[i] != RationalNumber(-1):
                        right_string += str(abs(self.coefficient[i]))
                        if not self.coefficient[i].isinteger():
                            right_string += " "
                    right_string += self.things[i].string + " + "
                else:
                    pass
            return left_string[:-3] + " == " + right_string[:-3]


try:
    a = input()
    equation_in = Equation(a)
    equation_in.solve()
    print(equation_in)
except:
    print("出现错误，请更正后重试")

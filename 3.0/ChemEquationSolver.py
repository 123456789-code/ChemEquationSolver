import copy


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

    def __repr__(self):
        return str(self)

    def isinteger(self):
        return self.q == 1


class Matrix:
    def __init__(self, lst, outcome):
        self.list = lst
        self.outcome = outcome
        self.upper_triangular_form_ed = False
        self.del_zero_row_ed = False
        self.reduced_row_echelon_form_ed = False
        self.major = []

    def __str__(self):
        string = ""
        for i in range(len(self.list)):
            for j in self.list[i]:
                string += str(j) + ", "
            string += "|| " + str(self.outcome[i]) + "\n"
        return string

    def __repr__(self):
        return str(self)

    def del_zero_row(self):
        lst = []
        for i in range(len(self.list)):
            if self.list[i] == len(self.list[i]) * [RationalNumber(0)] and self.outcome[i] == RationalNumber(0):
                lst.append(i)
        for i in range(len(lst)):
            lst[i] -= i
        for i in lst:
            del self.list[i]
            del self.outcome[i]
        return self

    def upper_triangular_form(self):  # 行阶梯矩阵
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
        self.del_zero_row()
        self.upper_triangular_form_ed = True
        return self

    def diagonal_eliminate(self):
        if not self.upper_triangular_form_ed:
            self.upper_triangular_form()
        if len(self.list[0]) != len(self.list):
            raise RuntimeError("this matrix can't be diagonal eliminated")
        for i in range(len(self.outcome)):
            if self.list[-1 - i][-1 - i] == RationalNumber(0):
                raise RuntimeError("there's a 0 on the diagonal")
            for j in range(len(self.outcome) - i - 1):
                factor = self.list[j][-1 - i] / self.list[-1 - i][-1 - i]
                for k in range(i, len(self.list[0])):
                    self.list[j][k] -= factor * self.list[-1 - i][k]
                self.outcome[j] -= factor * self.outcome[-1 - i]
        return self

    def solve(self):
        self.upper_triangular_form()
        # print(self)
        if len(self.list[0]) > len(self.list):
            return ["Infinity_solution"]
        elif len(self.list[0]) < len(self.list):
            return ["No_solution"]
        else:
            try:
                self.diagonal_eliminate()
            except RuntimeError("there's a 0 on the diagonal"):
                return ["Infinity_solution"]
            coefficient = []
            for i in range(len(self.outcome)):
                if self.list[i][i] == RationalNumber(0):
                    if self.outcome[i] == RationalNumber(0):
                        coefficient.append(RationalNumber(0))
                    else:
                        return ["No_solution"]
                else:
                    coefficient.append(self.outcome[i] / self.list[i][i])
            return coefficient

    def reduced_row_echelon_form(self):  # 会损坏outcome，请保证outcome全为0
        i, main_position = 0, 0
        while i < len(self.list) and main_position < len(self.list[0]):
            if self.list[i][main_position] == RationalNumber(0):
                for j in range(i + 1, len(self.list)):
                    if self.list[j][main_position] != RationalNumber(0):
                        self.list[i], self.list[j] = self.list[j], self.list[i]
                        break
                if self.list[i][main_position] == RationalNumber(0):
                    main_position += 1
                    continue
            for j in range(i + 1, len(self.list)):
                factor = self.list[j][main_position] / self.list[i][main_position]
                for k in range(main_position, len(self.list[0])):
                    self.list[j][k] -= factor * self.list[i][k]
            self.major.append(main_position)
            i += 1
            main_position += 1
        self.del_zero_row()
        for i in range(len(self.major)):
            if self.list[i][self.major[i]] != RationalNumber(1):
                for j in range(self.major[i] + 1, len(self.list[0])):
                    self.list[i][j] /= self.list[i][self.major[i]]
                self.list[i][self.major[i]] = RationalNumber(1)
            for j in range(i):
                factor = self.list[j][self.major[i]]
                for k in range(self.major[i], len(self.list[0])):
                    self.list[j][k] -= factor * self.list[i][k]
        self.reduced_row_echelon_form_ed = True
        return self

    def NullSpace(self):
        print(self)
        self.reduced_row_echelon_form()
        print(self)
        null_space = []
        free = list(set(range(len(self.list[0]))) - set(self.major))
        for i in range(len(free)):
            list_cache = copy.deepcopy(self.list)
            outcome_cache = copy.deepcopy(self.outcome)
            for j in range(len(free)):
                list_cache.append([RationalNumber(0)] * free[j] + [RationalNumber(1)] + [RationalNumber(0)] * (
                            len(self.list[0]) - free[j] - 1))
                if j == i:
                    outcome_cache.append(RationalNumber(-1))
                else:
                    outcome_cache.append(RationalNumber(0))
            matrix_cache = Matrix(list_cache, outcome_cache)
            null_space.append(matrix_cache.solve())
        return null_space


class Thing:
    def add_dict(self, d1, d2):
        for k, v in d2.items():
            if k in d1:
                d1[k] += v
            else:
                d1[k] = v
        return d1

    def times_dict(self, d, n):
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

    def __str__(self):
        return str(self.elements) + "\nElectron: " + self.electron

    def __repr__(self):
        return str(self)

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
                position += 1
                while position < len(string):
                    char = string[position]
                    if char.isdigit():
                        num_cache += char
                        position += 1
                    else:
                        element_dict = self.add_dict(element_dict,
                                                     self.times_dict(dict_cache,
                                                                     (1 if num_cache == "" else int(num_cache))))
                        dict_cache = {}
                        name_cache = ""
                        num_cache = ""
                        position -= 1
                        break
            elif char == ")":
                break
            elif char.isdigit():
                num_cache += char
            elif char.islower():
                name_cache += char
            elif char.isupper:
                element_dict = self.add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
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
        if "(" in element_dict:
            del element_dict["("]
        if ")" in element_dict:
            del element_dict[")"]
        print(element_dict)
        return element_dict, position


class Equation:
    def __init__(self, string):
        self.things = []
        self.coefficient_in = []
        self.coefficient = []
        self.flag = True  # True为未指定系数
        lst = string.split("--")
        for i in lst[0].split(" + "):
            i = i.replace(" ", "").replace("\n", "")
            if i[0].isdigit() or i[0] == "-":
                self.flag = False
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
                self.flag = False
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
        if self.flag:
            matrix1 = Matrix(matrix_cache, outcome_cache).del_zero_row()
            null_space = matrix1.NullSpace()
            # print(null_space)
            for i in range(len(null_space)):
                lcm = 1
                for j in null_space[i]:
                    lcm = lcm * j.q // j.gcd(lcm, j.q)
                lcm_ = RationalNumber(lcm)
                for j in range(len(null_space[i])):
                    null_space[i][j] *= lcm_
            self.coefficient = null_space
        else:
            for i in range(len(self.coefficient_in)):
                if self.coefficient_in[i] != RationalNumber(0):
                    matrix_cache.append([RationalNumber(0)] * i + [RationalNumber(1)] +
                                        [RationalNumber(0)] * (len(matrix_cache[0]) - i - 1))
                    outcome_cache.append(self.coefficient_in[i])
            matrix1 = Matrix(matrix_cache, outcome_cache).del_zero_row()
            self.coefficient = matrix1.solve()
            for i in range(len(self.coefficient)):
                if self.coefficient_in[i] != RationalNumber(0) and self.coefficient_in[i] != self.coefficient[i]:
                    self.coefficient = ["No_solution"]
                    break

    def single_equation_str(self, coefficient):
        left_string, right_string = "", ""
        for i in range(len(coefficient)):
            if coefficient[i] > RationalNumber(0):
                if coefficient[i] != RationalNumber(1):
                    left_string += str(coefficient[i])
                    if not coefficient[i].isinteger():
                        left_string += " "
                left_string += self.things[i].string + " + "
            elif coefficient[i] < RationalNumber(0):
                if coefficient[i] != RationalNumber(-1):
                    right_string += str(abs(coefficient[i]))
                    if not coefficient[i].isinteger():
                        right_string += " "
                right_string += self.things[i].string + " + "
            else:
                pass
        return left_string[:-3] + " == " + right_string[:-3]

    def __str__(self):
        if self.coefficient == ["No_solution"]:
            return "无解"
        elif self.coefficient == ["Infinity_solution"]:
            return "无穷解"
        else:
            if self.flag:
                string = ""
                for coefficient in self.coefficient:
                    string += self.single_equation_str(coefficient) + "\n"
                return string
            else:
                return self.single_equation_str(self.coefficient)

    def __repr__(self):
        return str(self)


# try:
a = input()
equation_in = Equation(a)
equation_in.solve()
print(equation_in)
# except:
#     print("出现错误，请更正后重试")

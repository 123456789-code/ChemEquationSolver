class RationalNumber:
    def __init__(self, p=0, q=1):
        self.p, self.q = p, q

    @staticmethod
    def gcd(a: int, b: int) -> int:
        if a < b:
            a, b = b, a
        if b == 0:
            return a if a != 0 else 1
        else:
            while a % b != 0:
                a, b = b, a % b
            return b

    def simplify(self):
        if self.q < 0:
            self.p, self.q = -self.p, -self.q
        gcd = self.gcd(abs(self.p), self.q)
        self.p //= gcd
        self.q //= gcd
        return self

    def __add__(self, other):
        return RationalNumber(self.p * other.q + self.q * other.p, self.q * other.q).simplify()

    def __sub__(self, other):
        return RationalNumber(self.p * other.q - self.q * other.p, self.q * other.q).simplify()

    def __mul__(self, other):
        return RationalNumber(self.p * other.p, self.q * other.q).simplify()

    def __truediv__(self, other):
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
        if self.q == 1:
            return self.p
        else:
            raise RuntimeError("rational_number in int() is not an integer")

    def __repr__(self):
        return str(self)

    def isinteger(self):
        return self.q == 1


R0 = RationalNumber(0)
R1 = RationalNumber(1)
RN1 = RationalNumber(-1)


class Matrix:
    def __init__(self, lst, outcome):
        self.list = lst
        self.outcome = outcome
        self.row_len = len(lst)
        self.column_len = len(lst[0])
        self.major = []

    def __str__(self):
        string = ""
        for i in range(self.row_len):
            for j in self.list[i]:
                string += str(j) + ", "
            string += "|| " + str(self.outcome[i]) + "\n"
        return string

    def __repr__(self):
        return str(self)

    @staticmethod
    def all_zero(lst):
        flag = True
        for i in lst:
            if i != R0:
                flag = False
                break
        return flag

    def count_list(self):
        self.row_len = len(self.list)
        self.column_len = len(self.list[0])

    def del_zero_row(self):
        while True:
            if self.all_zero(self.list[-1]) and self.outcome[-1] == R0:
                del self.list[-1]
                del self.outcome[-1]
            else:
                break
        self.count_list()

    def diagonal_form(self):
        for i in range(self.column_len):
            if i >= self.row_len:
                """Index out of range."""
                raise RuntimeError
            if self.list[i][i] == R0:
                for j in range(i + 1, self.row_len):
                    if self.list[j][i] != R0:
                        self.list[i], self.list[j] = self.list[j], self.list[i]
                        self.outcome[i], self.outcome[j] = self.outcome[j], self.outcome[i]
                        break
                if self.list[i][i] == R0:
                    """There's a 0 on the diagonal"""
                    raise RuntimeError
            for j in list(range(i)) + list(range(i + 1, self.row_len)):
                factor = self.list[j][i] / self.list[i][i]
                for k in range(i, self.column_len):
                    self.list[j][k] -= factor * self.list[i][k]
                self.outcome[j] -= factor * self.outcome[i]
        self.del_zero_row()
        return self

    def solve(self) -> [RationalNumber]:
        try:
            self.diagonal_form()
        except RuntimeError:
            return ["Infinity_solution"]
        coefficient = []
        for i in range(self.row_len):
            if self.list[i][i] == R0:
                if self.outcome[i] == R0:
                    return ["Infinity_solution"]
                else:
                    return ["No_solution"]
            else:
                coefficient.append(self.outcome[i] / self.list[i][i])
        return coefficient

    def reduced_row_echelon_form(self):
        if not self.all_zero(self.outcome):
            raise RuntimeError("row_echelon_form() can't be used when outcome includes not-zero rational number")
        i, main_position = 0, 0
        while i < self.row_len and main_position < self.column_len:
            if self.list[i][main_position] == R0:
                for j in range(i + 1, self.row_len):
                    if self.list[j][main_position] != R0:
                        self.list[i], self.list[j] = self.list[j], self.list[i]
                        break
                if self.list[i][main_position] == R0:
                    main_position += 1
                    continue
            for j in list(range(i)) + list(range(i + 1, self.row_len)):
                factor = self.list[j][main_position] / self.list[i][main_position]
                for k in range(main_position, self.column_len):
                    self.list[j][k] -= factor * self.list[i][k]
            if self.list[i][main_position] != R1:
                for k in range(main_position + 1, self.column_len):
                    self.list[i][k] /= self.list[i][main_position]
                self.list[i][main_position] = R1
            self.major.append(main_position)
            i += 1
            main_position += 1
        self.del_zero_row()
        return self

    def null_space(self) -> [[RationalNumber]]:
        self.reduced_row_echelon_form()
        null_space_return = []
        free = list(set(range(self.column_len)) - set(self.major))
        for i in range(len(free)):
            list_cache = []
            for j in self.list:
                list_cache.append(j[free[i]])
            list_cache += [R0] * i + [RN1] + [R0] * (len(free) - i - 1)
            null_space_return.append(list_cache)
        return null_space_return


class Thing:
    @staticmethod
    def add_dict(d1: dict, d2: dict) -> dict:
        for k, v in d2.items():
            if k in d1:
                d1[k] += v
            else:
                d1[k] = v
        return d1

    @staticmethod
    def times_dict(d: dict, n: int) -> dict:
        for i in d.keys():
            d[i] *= n
        return d

    def __init__(self, string):
        self.string = string.replace(" ", "").replace("\n", "")
        if self.string == "e[-]":
            self.elements = {"e": 0}
            self.electron = -1
        elif "]" in self.string:
            lst = self.string[:-1].split("[")
            if lst[1] == "":
                self.electron = 0
            elif lst[1] == "+":
                self.electron = 1
            elif lst[1] == "-":
                self.electron = -1
            else:
                self.electron = int(lst[1][-1] + lst[1][:-1])
            self.elements = self.separate(0)[0]
        else:
            self.electron = 0
            self.elements = self.separate(0)[0]

    def __str__(self):
        return str(self.elements) + "\nElectron: " + str(self.electron)

    def __repr__(self):
        return str(self)

    def separate(self, position):
        element_dict = {}
        name_cache = ""
        num_cache = ""
        dict_cache = {}
        while position < len(self.string):
            char = self.string[position]
            if char == "[":
                break
            elif char == "(":
                if name_cache != "":
                    element_dict = self.add_dict(element_dict, {name_cache: (1 if num_cache == "" else int(num_cache))})
                    name_cache = ""
                    num_cache = ""
                dict_cache, position = self.separate(position + 1)
                position += 1
                while position < len(self.string):
                    char = self.string[position]
                    if char.isdigit():
                        num_cache += char
                        position += 1
                    else:
                        element_dict = self.add_dict(
                            element_dict, self.times_dict(dict_cache, (1 if num_cache == "" else int(num_cache)))
                        )
                        dict_cache = {}
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
        element_dict = self.add_dict(
            element_dict, self.times_dict(dict_cache, (1 if num_cache == "" else int(num_cache)))
        )
        if "" in element_dict:
            del element_dict[""]
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
                self.coefficient_in.append(RationalNumber(i[:cache]))
                self.things.append(Thing(i[cache:]))
            else:
                self.coefficient_in.append(R0)
                self.things.append(Thing(i))
        for i in lst[1].split(" + "):
            i = i.replace(" ", "").replace("\n", "")
            if i[0].isdigit() or i[0] == "-":
                self.flag = False
                cache = 0
                while i[cache].isdigit() or i[cache] in "/-":
                    cache += 1
                self.coefficient_in.append(RationalNumber(i[:cache]) * RN1)
                self.things.append(Thing(i[cache:]))
            else:
                self.coefficient_in.append(R0)
                self.things.append(Thing(i))
        self.element_list = []
        for i in self.things:
            self.element_list += i.elements.keys()
        self.element_list = list(set(self.element_list))

    def solve(self):
        matrix_cache = list(list(RationalNumber(j.elements[i]) if i in j.elements else R0
                                 for j in self.things) for i in self.element_list)
        matrix_cache.append(list(RationalNumber(i.electron) for i in self.things))
        outcome_cache = [R0] * len(matrix_cache)
        if self.flag:
            matrix1 = Matrix(matrix_cache, outcome_cache)
            null_space = matrix1.null_space()
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
                if self.coefficient_in[i] != R0:
                    matrix_cache.append([R0] * i + [R1] + [R0] * (len(matrix_cache[0]) - i - 1))
                    outcome_cache.append(self.coefficient_in[i])
            matrix1 = Matrix(matrix_cache, outcome_cache)
            self.coefficient = matrix1.solve()
            for i in range(len(self.coefficient)):
                if self.coefficient_in[i] != self.coefficient[i]:
                    self.coefficient = ["No_solution"]
                    break

    def single_equation_str(self, coefficient):
        left_string, right_string = "", ""
        for i in range(len(coefficient)):
            if coefficient[i] > R0:
                if coefficient[i] != R1:
                    left_string += str(coefficient[i])
                    if not coefficient[i].isinteger():
                        left_string += " "
                left_string += self.things[i].string + " + "
            elif coefficient[i] < R0:
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
                return string.strip()
            else:
                return self.single_equation_str(self.coefficient)

    def __repr__(self):
        return str(self)


# try:
input_ = input()
equation_in = Equation(input_)
equation_in.solve()
print(equation_in)
# except:
#     print("出现错误，请更正后重试")

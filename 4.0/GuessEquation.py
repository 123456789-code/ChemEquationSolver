import os


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
                    element_dict = self.add_dict(element_dict, {name_cache: (int(num_cache) if num_cache else 1)})
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
                            element_dict, self.times_dict(dict_cache, (int(num_cache) if num_cache else 1))
                        )
                        dict_cache = {}
                        num_cache = ""
                        position -= 1
                        break
            elif char == ")":
                raise RuntimeError("括号不匹配")
            elif char.isdigit():
                num_cache += char
            elif char.islower():
                name_cache += char
            elif char.isupper:
                element_dict = self.add_dict(element_dict, {name_cache: (int(num_cache) if num_cache else 1)})
                name_cache = char
                num_cache = ""
            else:
                raise RuntimeError("非法字符")
            position += 1
        element_dict = self.add_dict(element_dict, {name_cache: (int(num_cache) if num_cache else 1)})
        element_dict = self.add_dict(
            element_dict, self.times_dict(dict_cache, (int(num_cache) if num_cache else 1))
        )
        if "" in element_dict:
            del element_dict[""]
        return element_dict, position


def guess_equation(lst):
    if len(lst) == 1:
        return "一个物质反应个屁啊"
    elif len(lst) == 0:
        return "你tm一个物质都没有"
    else:
        things = []
        for i in lst:
            things.append(Thing(i))
        element_list = []
        for i in things:
            element_list += i.elements.keys()
        element_list = list(set(element_list))
        matrix_cache = list(list(RationalNumber(j.elements[i]) if i in j.elements else R0
                                 for j in things) for i in element_list)
        matrix_cache.append(list(RationalNumber(i.electron) for i in things))
        outcome_cache = [R0] * len(matrix_cache)
        matrix1 = Matrix(matrix_cache, outcome_cache)
        null_space = matrix1.null_space()
        for i in range(len(null_space)):
            lcm = 1
            for j in null_space[i]:
                lcm = lcm * j.q // j.gcd(lcm, j.q)
            lcm_ = RationalNumber(lcm)
            for j in range(len(null_space[i])):
                null_space[i][j] *= lcm_
        if not null_space:
            return "无反应发生"
        string = ""
        for coefficient in null_space:
            left_string, right_string = "", ""
            for i in range(len(coefficient)):
                if coefficient[i] > R0:
                    if coefficient[i] != R1:
                        left_string += str(coefficient[i])
                        if not coefficient[i].isinteger():
                            left_string += " "
                    left_string += things[i].string + " + "
                elif coefficient[i] < R0:
                    if coefficient[i] != RN1:
                        right_string += str(abs(coefficient[i]))
                        if not coefficient[i].isinteger():
                            right_string += " "
                    right_string += things[i].string + " + "
                else:
                    pass
            string += left_string[:-3] + " == " + right_string[:-3] + "\n"
        string = string.strip()
        return string


if __name__ == "__main__":
    print("请输入物质:")
    lst = []
    string = "none"
    while string != "":
        string = input().strip()
        lst.append(string)
    lst.pop()
    print(guess_equation(lst))
    os.system("pause>nul")

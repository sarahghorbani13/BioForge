from datetime import datetime


# کلاس اصلی برای کار روی توالی‌های دی‌ ان‌ ای
# ورودی: دیکشنری رکوردها (کلید = شناسه، مقدار = توالی) که پارسر فایل فستا ساخته
class DNASequence:

    # سازنده: توالی‌ها را از دیکشنری بیرون می‌کشد و بزرگ‌ حرف می‌کند
    def __init__(self, sequence):
        # فقط مقدارها (خود توالی‌ها) را نگه می‌داریم و همه را بزرگ‌ حرف می‌کنیم
        self.sequence = [seq.upper() for seq in list(sequence.values())]
        self.valid_bases = {"A", "T", "C", "G"}
        # بلافاصله بعد از ساخت، توالی‌ها را بررسی می‌کنیم
        self.validation()

    # اعتبارسنجی: اگر حرف غیرمجاز باشد، خطا را در فایل لاگ ثبت می‌کند
    def validation(self):
        try:
            for dna in self.sequence:
                # اگر مجموعه‌ی حروف توالی زیرمجموعه‌ی بازهای مجاز نباشد، یعنی حرف اشتباه دارد
                if not set(dna).issubset(self.valid_bases):
                    raise Exception("InvalidSequenceError")
        except Exception as dna_operation_error:
            # خطا را همراه زمان در فایل لاگ اضافه می‌کنیم (حالت افزودن یعنی به انتهای فایل می‌نویسد)
            with open("bioforge.log", "a") as file:
                file.write(f"{datetime.now()} - {dna_operation_error}\n")

    #مکمل 
    def complement(self):
        # از لیست اصلی کپی می‌گیریم تا توالی‌های اصلی خراب نشوند
        complemented = self.sequence.copy()

        for i in range(len(complemented)):
            # رشته را به لیست حروف تبدیل می‌کنیم تا بشود تک‌ تک تغییر داد
            dna = list(complemented[i])

            for j in range(len(dna)):
                if dna[j] == 'A':
                    dna[j] = 'T'
                elif dna[j] == 'T':
                    dna[j] = 'A'
                elif dna[j] == 'C':
                    dna[j] = 'G'
                elif dna[j] == 'G':
                    dna[j] = 'C'

            # حروف را دوباره به یک رشته تبدیل می‌کنیم
            complemented[i] = "".join(dna)

        return complemented

    # معکوس‌ مکمل: اول مکمل می‌گیریم، بعد ترتیب را برعکس می‌کنیم
    def reverse_complement(self):
        rev_complemented = self.complement()

        for i in range(len(rev_complemented)):
            # برش با گام منفی یک یعنی رشته را از آخر به اول بخوان
            rev_complemented[i] = rev_complemented[i][::-1]

        return rev_complemented

    # تبدیل دی‌ ان‌ ای به آر‌ ان‌ ای
    def dna_to_rna(self):
        d_to_r = self.sequence.copy()

        for i in range(len(d_to_r)):
            d_to_r[i] = d_to_r[i].replace("T", "U")

        return d_to_r

    # درصد جی‌ سی:(میانگین روی همه‌ی توالی‌ها)
    def gc_content(self):
        total_gc_percent = 0

        for seq in self.sequence:
            c_count = seq.count("C")
            g_count = seq.count("G")
            # درصد سیتوزین و گوانین برای همین توالی
            total_gc_percent += ((c_count + g_count) / len(seq)) * 100

        # میانگین درصدها روی تعداد توالی‌ها
        return total_gc_percent / len(self.sequence)

    # خروجی نهایی این کلاس برای مرحله‌ی بعد
    # برای هر توالی دو رشته می‌دهد: مستقیم و معکوس‌مکمل
    # هر مورد یک جفت (رشته‌ی آر‌ان‌ای، جهت رشته) است
    def rna_strands(self):
        rna_strands = []

        # رشته‌ی مستقیم: همان آر‌ان‌ای هر توالی
        for rna in self.dna_to_rna():
            rna_strands.append((rna, "forward"))

        # رشته‌ی ریورس کامپلمنت: معکوس‌ مکمل را می‌گیریم و تی را به یو تبدیل می‌کنیم
        for dna in self.reverse_complement():
            rna = dna.replace("T", "U")
            rna_strands.append((rna, "reverse_complement"))

        return rna_strands

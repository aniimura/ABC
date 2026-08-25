import ABC3.Found.GaloisRep.TateLeading

/-!
# Galois (G6) 第 259 ブロック —— **★★★★★★★★★葉 (d) が終わった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> `a·w = q = a'·w'`、`a, w, a', w' ∈ I` のとき
> `X(a,w) = X(a',w')` かつ `Y(a,w) = Y(a',w')` なら `(a,w) = (a',w')`(`tate_inj`)。

さらに `R` が整域なら **`X` だけで `±` を除いて決まる**(`tate_inj_X`):

    X(a,w) = X(a',w')  ⟹  (a,w) = (a',w')  または  (a,w) = (w',a')

★後者は `−P` に対応する(第 213 の `tateYpair_eq_negY`)。

## ★★★★★★★逐次近似 —— 付値を使わない

    a − a' ∈ I^j  かつ  w − w' ∈ I^j   ⟹   a − a' ∈ I^{j+1} かつ w − w' ∈ I^{j+1}

を `j` について繰り返し、`IsHausdorff`(`eq_zero_of_mem_pow`)で `a = a'`。
★★**付値も基本領域も使わない**——`I` 進の言葉だけで閉じる。

### ★★★1 段の中身

`q` は両側で同じなので、`X` の差でも `Y` の差でも **`s₁(q)` の項が消える**。残るのは

| 差 | 位 |
|---|---|
| `g(a) − g(a')`、`g(w) − g(w')` | `I^{j+1}`(`g` は 2 次から始まる) |
| `f(w) − f(w') − (w − w')` | `I^{j+1}`(第 258) |
| 尾の差 | `I^{j+2}`(`q ∈ I²` を 1 つ使う) |

★`Y` の差から **`w − w' ∈ I^{j+1}`** が出る(`Y ≡ −w`)。
それを `X` の差(`X ≡ a + w`)に入れると **`a − a' ∈ I^{j+1}`** が出る。
★★`Y` を先に、`X` を後に、という順序が要点である。

## ★★`X` だけの版は曲線の式から

`X = X'` なら曲線の式(葉 (b)、第 240)の差から

    (Y − Y')·(Y + Y' + X) = 0

`R` が整域なので分岐は 2 つだけ。第 2 の枝は `Y' = −X − Y`、すなわち
`(w', a')` が `(a, w)` と同じ座標を与える(`tateXpair_symm`・`tateYpair_swap`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `adicSum_diff_mem` | ★★★★adic 和の差 |
| `tateYterm_diff_mem` | ★★★★★`g` の差 |
| `tateXtail_diff_mem`・`tateYtail_diff_mem` | ★★★★尾の差 |
| `tate_inj_step` | ★★★★★★★**単射性の 1 段** |
| `tate_inj` | ★★★★★★★★★**葉 (d)** |
| `tate_inj_X` | ★★★★★★★★`X` だけの版 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★差の道具 -/

/-- ★★★★**adic 和の差**——各項の差が `I^k` なら和の差も `I^k`。 -/
theorem adicSum_diff_mem [IsPrecomplete I R] (a b : ℕ → R) (ha : ∀ n, a n ∈ I ^ n)
    (hb : ∀ n, b n ∈ I ^ n) (k : ℕ) (hk : ∀ n, a n - b n ∈ I ^ k) :
    adicSum a ha - adicSum b hb ∈ I ^ k := by
  have hA := (smodEq_iff_sub_mem (I := I) (n := k) (partialSum a k) (adicSum a ha)).1
    (adicSum_spec a ha k)
  have hB := (smodEq_iff_sub_mem (I := I) (n := k) (partialSum b k) (adicSum b hb)).1
    (adicSum_spec b hb k)
  have hpart : partialSum a k - partialSum b k ∈ I ^ k := by
    rw [partialSum, partialSum, ← Finset.sum_sub_distrib]
    exact Ideal.sum_mem _ fun n _ => hk n
  have hkey : adicSum a ha - adicSum b hb
      = -(partialSum a k - adicSum a ha) + (partialSum b k - adicSum b hb)
        + (partialSum a k - partialSum b k) := by ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (neg_mem hA) hB) hpart

/-- ★★★★★**`g` の差**——`g(a) − g(b)` は 1 つ位が上がる。

★`D := (1−a)³(1−b)³` として `D·(g(a) − g(b)) = (a−b)(a + b − 3ab + a²b²)` である。 -/
theorem tateYterm_diff_mem [IsAdicComplete I R] {j : ℕ} {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (hab : a - b ∈ I ^ j) : tateYterm a - tateYterm b ∈ I ^ (j + 1) := by
  have hua : IsUnit (1 - a) := isUnit_one_sub ha
  have hub : IsUnit (1 - b) := isUnit_one_sub hb
  set D : R := (1 - a) ^ 3 * (1 - b) ^ 3 with hDdef
  have huD : IsUnit D := by
    simp only [hDdef]
    exact (hua.pow 3).mul (hub.pow 3)
  have hda : (1 - a) ^ 3 * tateYterm a = a ^ 2 := mul_tateYterm ha
  have hdb : (1 - b) ^ 3 * tateYterm b = b ^ 2 := mul_tateYterm hb
  have hD : D * (tateYterm a - tateYterm b)
      = (a - b) * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2) := by
    simp only [hDdef]
    linear_combination (1 - b) ^ 3 * hda - (1 - a) ^ 3 * hdb
  have hP : a + b - 3 * (a * b) + a ^ 2 * b ^ 2 ∈ I := by
    refine Ideal.add_mem _ (Ideal.sub_mem _ (Ideal.add_mem _ ha hb)
      (Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ ha))) ?_
    exact Ideal.mul_mem_right _ _ (by rw [sq]; exact Ideal.mul_mem_right _ _ ha)
  have h1 : Ring.inverse D * D = 1 := Ring.inverse_mul_cancel _ huD
  have hkey : tateYterm a - tateYterm b
      = (a - b) * (Ring.inverse D * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by
    calc tateYterm a - tateYterm b
        = Ring.inverse D * (D * (tateYterm a - tateYterm b)) := by
          rw [← mul_assoc, h1, one_mul]
      _ = Ring.inverse D * ((a - b) * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by rw [hD]
      _ = (a - b) * (Ring.inverse D * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by ring
  rw [hkey, pow_succ]
  exact Ideal.mul_mem_mul hab (Ideal.mul_mem_left _ _ hP)

theorem tateXterm_diff_mem' [IsAdicComplete I R] {j : ℕ} {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (hab : a - b ∈ I ^ j) : tateXterm a - tateXterm b ∈ I ^ j := by
  have h := tateXterm_diff_mem ha hb hab
  have h2 := Ideal.add_mem (I ^ j) (Ideal.pow_le_pow_right (Nat.le_succ j) h) hab
  simpa using h2

theorem tateYterm_diff_mem' [IsAdicComplete I R] {j : ℕ} {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (hab : a - b ∈ I ^ j) : tateYterm a - tateYterm b ∈ I ^ j :=
  Ideal.pow_le_pow_right (Nat.le_succ j) (tateYterm_diff_mem ha hb hab)

theorem pow_succ_mul_mem' {q : R} (hq : q ∈ I) (n : ℕ) (u : R) : q ^ (n + 1) * u ∈ I := by
  rw [pow_succ]
  exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ hq)

theorem mul_pow_diff_mem {j m n : ℕ} {a b q : R} (hqm : q ∈ I ^ m) (hab : a - b ∈ I ^ j) :
    q ^ (n + 1) * a - q ^ (n + 1) * b ∈ I ^ (j + m) := by
  have h1 : (a - b) * q ∈ I ^ (j + m) := by
    rw [pow_add]
    exact Ideal.mul_mem_mul hab hqm
  have h2 : q ^ (n + 1) * a - q ^ (n + 1) * b = q ^ n * ((a - b) * q) := by ring
  rw [h2]
  exact Ideal.mul_mem_left _ _ h1

theorem tateXtail_diff_mem [IsAdicComplete I R] {j m : ℕ} (a b q : R) (hq : q ∈ I)
    (hqm : q ∈ I ^ m) (hab : a - b ∈ I ^ j) :
    tateXtail a q hq - tateXtail b q hq ∈ I ^ (j + m) := by
  refine adicSum_diff_mem _ _ _ _ (j + m) fun n => ?_
  exact tateXterm_diff_mem' (pow_succ_mul_mem' hq n a) (pow_succ_mul_mem' hq n b)
    (mul_pow_diff_mem hqm hab)

theorem tateYtail_diff_mem [IsAdicComplete I R] {j m : ℕ} (a b q : R) (hq : q ∈ I)
    (hqm : q ∈ I ^ m) (hab : a - b ∈ I ^ j) :
    tateYtail a q hq - tateYtail b q hq ∈ I ^ (j + m) := by
  refine adicSum_diff_mem _ _ _ _ (j + m) fun n => ?_
  exact tateYterm_diff_mem' (pow_succ_mul_mem' hq n a) (pow_succ_mul_mem' hq n b)
    (mul_pow_diff_mem hqm hab)

/-! ## ★★★★★★★逐次近似 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**単射性の 1 段**——`I^j` の一致は `I^{j+1}` に上がる。

★`Y` の差から `w − w'` が、それを使って `X` の差から `a − a'` が出る(順序が要点)。 -/
theorem tate_inj_step [IsAdicComplete I R] {j : ℕ} (a w a' w' q : R) (hq : q ∈ I)
    (ha : a ∈ I) (hw : w ∈ I) (ha' : a' ∈ I) (hw' : w' ∈ I)
    (haw : a * w = q)
    (hX : tateXpair a w q hq = tateXpair a' w' q hq)
    (hY : tateYpair a w q hq = tateYpair a' w' q hq)
    (hja : a - a' ∈ I ^ j) (hjw : w - w' ∈ I ^ j) :
    a - a' ∈ I ^ (j + 1) ∧ w - w' ∈ I ^ (j + 1) := by
  have hq2 : q ∈ I ^ 2 := by
    rw [← haw, pow_two]
    exact Ideal.mul_mem_mul ha hw
  have hup : ∀ x : R, x ∈ I ^ (j + 2) → x ∈ I ^ (j + 1) := fun x hx =>
    Ideal.pow_le_pow_right (by omega) hx
  have dGA : tateYterm a - tateYterm a' ∈ I ^ (j + 1) := tateYterm_diff_mem ha ha' hja
  have dGW : tateYterm w - tateYterm w' ∈ I ^ (j + 1) := tateYterm_diff_mem hw hw' hjw
  have dFW : tateXterm w - tateXterm w' - (w - w') ∈ I ^ (j + 1) :=
    tateXterm_diff_mem hw hw' hjw
  have dFA : tateXterm a - tateXterm a' - (a - a') ∈ I ^ (j + 1) :=
    tateXterm_diff_mem ha ha' hja
  have dTGA : tateYtail a q hq - tateYtail a' q hq ∈ I ^ (j + 1) :=
    hup _ (tateYtail_diff_mem a a' q hq hq2 hja)
  have dTFA : tateXtail a q hq - tateXtail a' q hq ∈ I ^ (j + 1) :=
    hup _ (tateXtail_diff_mem a a' q hq hq2 hja)
  have dTFW : tateXtail w q hq - tateXtail w' q hq ∈ I ^ (j + 1) :=
    hup _ (tateXtail_diff_mem w w' q hq hq2 hjw)
  have dTGW : tateYtail w q hq - tateYtail w' q hq ∈ I ^ (j + 1) :=
    hup _ (tateYtail_diff_mem w w' q hq hq2 hjw)
  have hY0 : tateYpair a w q hq - tateYpair a' w' q hq = 0 := by rw [hY]; ring
  have hX0 : tateXpair a w q hq - tateXpair a' w' q hq = 0 := by rw [hX]; ring
  have hWkey : w - w' = (tateYterm a - tateYterm a') + (tateYtail a q hq - tateYtail a' q hq)
      - (tateXterm w - tateXterm w' - (w - w')) - (tateXtail w q hq - tateXtail w' q hq)
      - (tateYterm w - tateYterm w') - (tateYtail w q hq - tateYtail w' q hq)
      - (tateYpair a w q hq - tateYpair a' w' q hq) := by
    rw [tateYpair, tateYpair]
    ring
  have hWmem : w - w' ∈ I ^ (j + 1) := by
    rw [hWkey, hY0, sub_zero]
    exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
      (Ideal.add_mem _ dGA dTGA) dFW) dTFW) dGW) dTGW
  refine ⟨?_, hWmem⟩
  have hXkey : a - a' = (tateXpair a w q hq - tateXpair a' w' q hq)
      - (tateXterm a - tateXterm a' - (a - a')) - (tateXtail a q hq - tateXtail a' q hq)
      - (w - w') - (tateXterm w - tateXterm w' - (w - w'))
      - (tateXtail w q hq - tateXtail w' q hq) := by
    rw [tateXpair, tateXpair]
    ring
  rw [hXkey, hX0]
  exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (Ideal.zero_mem _) dFA) dTFA) hWmem) dFW) dTFW

/-! ## ★★★★★★★★★葉 (d) -/

/-- ★★★★★★★★★**葉 (d) —— Tate 座標は単射**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_inj [IsAdicComplete I R] (a w a' w' q : R) (hq : q ∈ I)
    (ha : a ∈ I) (hw : w ∈ I) (ha' : a' ∈ I) (hw' : w' ∈ I)
    (haw : a * w = q)
    (hX : tateXpair a w q hq = tateXpair a' w' q hq)
    (hY : tateYpair a w q hq = tateYpair a' w' q hq) :
    a = a' ∧ w = w' := by
  have hstep : ∀ j : ℕ, a - a' ∈ I ^ j ∧ w - w' ∈ I ^ j := by
    intro j
    induction j with
    | zero => simp
    | succ n ih => exact tate_inj_step a w a' w' q hq ha hw ha' hw' haw hX hY ih.1 ih.2
  refine ⟨?_, ?_⟩
  · exact sub_eq_zero.1 (eq_zero_of_mem_pow (I := I) fun j => (hstep j).1)
  · exact sub_eq_zero.1 (eq_zero_of_mem_pow (I := I) fun j => (hstep j).2)

/-- ★★★★★★★★**`X` だけでも `±` を除いて決まる**——`R` が整域のとき。

★曲線の式(葉 (b))の差から `(Y − Y')(Y + Y' + X) = 0`。 -/
theorem tate_inj_X [IsAdicComplete I R] [IsDomain R] (a w a' w' q : R) (hq : q ∈ I)
    (ha : a ∈ I) (hw : w ∈ I) (ha' : a' ∈ I) (hw' : w' ∈ I)
    (haw : a * w = q) (haw' : a' * w' = q)
    (hX : tateXpair a w q hq = tateXpair a' w' q hq) :
    (a = a' ∧ w = w') ∨ (a = w' ∧ w = a') := by
  have e1 := tate_equation a w q hq haw (isUnit_one_sub ha) (isUnit_one_sub hw)
  have e2 := tate_equation a' w' q hq haw' (isUnit_one_sub ha') (isUnit_one_sub hw')
  rw [← hX] at e2
  have hfac : (tateYpair a w q hq - tateYpair a' w' q hq)
      * (tateYpair a w q hq + tateYpair a' w' q hq + tateXpair a w q hq) = 0 := by
    linear_combination e1 - e2
  rcases mul_eq_zero.1 hfac with h | h
  · left
    exact tate_inj a w a' w' q hq ha hw ha' hw' haw hX (sub_eq_zero.1 h)
  · right
    have hX2 : tateXpair a w q hq = tateXpair w' a' q hq := by
      rw [tateXpair_symm a' w' q hq]; exact hX
    have hY2 : tateYpair a w q hq = tateYpair w' a' q hq := by
      rw [tateYpair_swap a' w' q hq, ← hX]
      linear_combination h
    exact tate_inj a w w' a' q hq ha hw hw' ha' haw hX2 hY2

/-! ## ★出典の紐付け(`.src`) -/

def tate_inj_step.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単射性の 1 段)",
    sectionId := "genell-def-3-3" }

def tate_inj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (d)、座標は単射)",
    sectionId := "genell-def-3-3" }

def tate_inj_X.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X だけで ± を除いて決まる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep

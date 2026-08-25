import ABC3.Found.GaloisRep.TatePointGroup

/-!
# Galois (G6) 第 262 ブロック —— **★★★★★★★葉 (e) の線型化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★`X + Y` が `a` を線型に取り出す

葉 (e)(全射性)は「点 `(x,y)` から `u` を作る」ことである。素朴に `X(u) = x` を
解こうとすると、`f(t) = t/(1−t)²` が **2 対 1**(`f(t) = f(1/t)`)なので
Hensel の判別式が出てきて残余体の条件に化ける。

★★★しかし `X` と `Y` を**両方**使うと線型になる:

    f(t) + g(t) = t/(1−t)³,     g(t) = t²/(1−t)³
    ⟹ g(t) = t·(f(t) + g(t))

すなわち **`Y = a·(X + Y)`**(主要項で)。`X + Y` は単元なので

    a = Y·(X + Y)⁻¹ + (補正)

★これが `a` を取り出す式である。**平方根も Hensel も要らない**。

## ★★誤差はすべて `I` の元

    Y − a(X + Y) = (Y − g(a)) − a·((X + Y) − (f(a) + g(a)))

★右辺の 2 つの括弧は、`w ∈ I` なら**尾と `s₁(q)` だけ**でできているので `I` の元である
(`tateYpair_sub_tateYterm_mem`・`tateXpair_add_tateYpair_sub_mem`)。
第 3 項 `g(a) − a(f(a)+g(a))` は**恒等的に 0**(`tateYterm_eq_mul`)。

## ★★単元 + `I` は単元

`IsAdicComplete I R` なら `1 − i` が単元(`isUnit_one_sub`)なので、
`u + i = u(1 − (−u⁻¹i))` も単元である(`isUnit_add_mem`)。
★これで `X + Y` の単元性が主要項 `a/(1−a)³` から出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mul_tateXYterm` | ★★★★★★`(1−t)³(f+g) = t` |
| `tateYterm_eq_mul` | ★★★★★**`g(a) = a(f(a)+g(a))`** |
| `tateYpair_sub_mul_mem` | ★★★★★★★**`Y − a(X+Y) ∈ I`** |
| `isUnit_add_mem` | ★★★★単元 + `I` は単元 |
| `isUnit_tateXpair_add_tateYpair` | ★★★★★★`X + Y` は単元 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★線型化の恒等式 -/

theorem mul_tateXterm' {t : R} (hu : IsUnit (1 - t)) : (1 - t) ^ 2 * tateXterm t = t := by
  rw [tateXterm,
    show (1 - t) ^ 2 * (t * Ring.inverse (1 - t) ^ 2)
      = t * ((1 - t) * Ring.inverse (1 - t)) ^ 2 by ring,
    Ring.mul_inverse_cancel _ hu, one_pow, mul_one]

theorem mul_tateYterm' {t : R} (hu : IsUnit (1 - t)) : (1 - t) ^ 3 * tateYterm t = t ^ 2 := by
  rw [tateYterm,
    show (1 - t) ^ 3 * (t ^ 2 * Ring.inverse (1 - t) ^ 3)
      = t ^ 2 * ((1 - t) * Ring.inverse (1 - t)) ^ 3 by ring,
    Ring.mul_inverse_cancel _ hu, one_pow, mul_one]

/-- ★★★★★★**`f(t) + g(t) = t/(1−t)³`**——`X + Y` は `t` について 1 次で始まる。 -/
theorem mul_tateXYterm {t : R} (hu : IsUnit (1 - t)) :
    (1 - t) ^ 3 * (tateXterm t + tateYterm t) = t := by
  have h1 := mul_tateXterm' hu
  have h2 := mul_tateYterm' hu
  have h3 : (1 - t) ^ 3 * tateXterm t = (1 - t) * t := by
    calc (1 - t) ^ 3 * tateXterm t = (1 - t) * ((1 - t) ^ 2 * tateXterm t) := by ring
      _ = (1 - t) * t := by rw [h1]
  calc (1 - t) ^ 3 * (tateXterm t + tateYterm t)
      = (1 - t) ^ 3 * tateXterm t + (1 - t) ^ 3 * tateYterm t := by ring
    _ = (1 - t) * t + t ^ 2 := by rw [h3, h2]
    _ = t := by ring

/-- ★★★★★**`g(a) = a·(f(a) + g(a))`**——`a` を線型に取り出す鍵。

★平方根も Hensel も要らないのはこの恒等式のおかげである。 -/
theorem tateYterm_eq_mul {a : R} (hu : IsUnit (1 - a)) :
    tateYterm a = a * (tateXterm a + tateYterm a) := by
  have h1 : Ring.inverse ((1 - a) ^ 3) * (1 - a) ^ 3 = 1 :=
    Ring.inverse_mul_cancel _ (hu.pow 3)
  have hsum := mul_tateXYterm hu
  have hg := mul_tateYterm' hu
  have hkey : (1 - a) ^ 3 * (tateYterm a - a * (tateXterm a + tateYterm a)) = 0 := by
    rw [mul_sub, hg, show (1 - a) ^ 3 * (a * (tateXterm a + tateYterm a))
      = a * ((1 - a) ^ 3 * (tateXterm a + tateYterm a)) by ring, hsum]
    ring
  have h2 := congrArg (fun z => Ring.inverse ((1 - a) ^ 3) * z) hkey
  simp only [← mul_assoc, h1, one_mul, mul_zero] at h2
  exact sub_eq_zero.1 h2

/-! ## ★★誤差はすべて `I` の元 -/

theorem tateXterm_mem_one {w : R} (hw : w ∈ I) : tateXterm w ∈ I := by
  have h := tateXterm_mem_pow (I := I) (k := 1) (by rwa [pow_one])
  rwa [pow_one] at h

theorem tateYterm_mem_one {w : R} (hw : w ∈ I) : tateYterm w ∈ I := by
  have h := tateYterm_mem_pow (I := I) (k := 1) (by rwa [pow_one])
  rwa [pow_one] at h

/-- ★★★★`Y` の主要項は `g(a)`。 -/
theorem tateYpair_sub_tateYterm_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (hw : w ∈ I) : tateYpair a w q hq - tateYterm a ∈ I := by
  have hkey : tateYpair a w q hq - tateYterm a
      = tateYtail a q hq + (-1) * tateXterm w + (-1) * tateXtail w q hq
        + (-1) * tateYterm w + (-1) * tateYtail w q hq
        + evalAdic (sigmaSeries 1) q hq := by
    rw [tateYpair]
    ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (tateYtail_mem a q hq) (Ideal.mul_mem_left _ _ (tateXterm_mem_one hw)))
    (Ideal.mul_mem_left _ _ (tateXtail_mem w q hq)))
    (Ideal.mul_mem_left _ _ (tateYterm_mem_one hw)))
    (Ideal.mul_mem_left _ _ (tateYtail_mem w q hq)))
    (evalAdic_mem (sigmaSeries 1) (by simp [sigmaSeries]) q hq)

/-- ★★★★`X + Y` の主要項は `f(a) + g(a)`。 -/
theorem tateXpair_add_tateYpair_sub_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (hw : w ∈ I) :
    tateXpair a w q hq + tateYpair a w q hq - (tateXterm a + tateYterm a) ∈ I := by
  have hkey : tateXpair a w q hq + tateYpair a w q hq - (tateXterm a + tateYterm a)
      = tateXtail a q hq + tateYtail a q hq + (-1) * tateYterm w + (-1) * tateYtail w q hq
        + (-1) * evalAdic (sigmaSeries 1) q hq := by
    rw [tateXpair, tateYpair]
    ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (Ideal.add_mem _ (tateXtail_mem a q hq) (tateYtail_mem a q hq))
    (Ideal.mul_mem_left _ _ (tateYterm_mem_one hw)))
    (Ideal.mul_mem_left _ _ (tateYtail_mem w q hq)))
    (Ideal.mul_mem_left _ _ (evalAdic_mem (sigmaSeries 1) (by simp [sigmaSeries]) q hq))

/-- ★★★★★★★**`a` は `Y/(X+Y)` で取り出せる**——差は `I` の元。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateYpair_sub_mul_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (hu : IsUnit (1 - a)) :
    tateYpair a w q hq - a * (tateXpair a w q hq + tateYpair a w q hq) ∈ I := by
  have h1 := tateYpair_sub_tateYterm_mem a w q hq hw
  have h2 := tateXpair_add_tateYpair_sub_mem a w q hq hw
  have h4 : tateYterm a - a * (tateXterm a + tateYterm a) = 0 :=
    sub_eq_zero.2 (tateYterm_eq_mul hu)
  have hkey : tateYpair a w q hq - a * (tateXpair a w q hq + tateYpair a w q hq)
      = (tateYpair a w q hq - tateYterm a)
        + (-a) * (tateXpair a w q hq + tateYpair a w q hq - (tateXterm a + tateYterm a))
        + (tateYterm a - a * (tateXterm a + tateYterm a)) := by ring
  rw [hkey, h4, add_zero]
  exact Ideal.add_mem _ h1 (Ideal.mul_mem_left _ _ h2)

/-! ## ★★単元性 -/

/-- ★★★★**単元 + `I` の元は単元**——`isUnit_one_sub` の言い換え。 -/
theorem isUnit_add_mem [IsAdicComplete I R] {u i : R} (hu : IsUnit u) (hi : i ∈ I) :
    IsUnit (u + i) := by
  obtain ⟨v, rfl⟩ := hu
  have hvi : (v : R) * ((v⁻¹ : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have h1 : ((v : R) + i) = (v : R) * (1 - (-(((v⁻¹ : Rˣ) : R) * i))) := by
    rw [mul_sub, mul_one, mul_neg, ← mul_assoc, hvi, one_mul]
    ring
  rw [h1]
  exact v.isUnit.mul (isUnit_one_sub (I := I) (neg_mem (Ideal.mul_mem_left _ _ hi)))

/-- ★★★★★★**`X + Y` は単元**——`a` が単元で `1 − a` が単元のとき。 -/
theorem isUnit_tateXpair_add_tateYpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (hw : w ∈ I) (ha : IsUnit a) (hu : IsUnit (1 - a)) :
    IsUnit (tateXpair a w q hq + tateYpair a w q hq) := by
  have hmul : IsUnit ((1 - a) ^ 3 * (tateXterm a + tateYterm a)) := by
    rw [mul_tateXYterm hu]; exact ha
  have hfg : IsUnit (tateXterm a + tateYterm a) := isUnit_of_mul_isUnit_right hmul
  have hd := tateXpair_add_tateYpair_sub_mem a w q hq hw
  have hsplit : tateXpair a w q hq + tateYpair a w q hq
      = (tateXterm a + tateYterm a)
        + (tateXpair a w q hq + tateYpair a w q hq - (tateXterm a + tateYterm a)) := by ring
  rw [hsplit]
  exact isUnit_add_mem hfg hd

/-! ## ★出典の紐付け(`.src`) -/

def tateYterm_eq_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——g(a) = a(f(a)+g(a)))",
    sectionId := "genell-def-3-3" }

def tateYpair_sub_mul_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——a は Y/(X+Y) で取り出せる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep

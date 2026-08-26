import ABC3.Found.GaloisRep.TateEquationMod

/-!
# Galois (G6) 第 213 ブロック —— **★★★★★★Tate 座標は `Kˣ/q^ℤ` 上で定まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★葉 (d) の前半——写像が商の上で定まること

Tate 一意化 `E(K) ≅ Kˣ/q^ℤ` の写像は `u` から作るが、**`u` と `qu` は同じ点に行く**
必要がある。★第 106 ブロックの基本領域(`0 ≤ v(u) < v(q)` の正規化)を使うと、
これは**代表元の一意性**に還元される。

### ★★★★代表元を関数にする

`exists_unique_normalized_rep`(第 106)は `∃!` だったので、`choose` で関数にした:

    normRep v q hq c : Kˣ        (c : Kˣ ⧸ q^ℤ)

★`mk (normRep c) = c`・`0 ≤ v(normRep c) < v(q)`、そして
**`eq_normRep`——同じ類の正規化元は代表元に一致する**。

### ★★★★★対 `(a, w)` は類だけで決まる

`R → K` が単射で `R` が整域なら:

| 段 | 内容 |
|---|---|
| 1 | 同じ類の正規化元は一致する(`eq_normRep`) |
| 2 | `algebraMap a = u = u' = algebraMap a'` と単射性で `a = a'` |
| 3 | `a·w = q = a'·w'` と `a ≠ 0` の相殺で `w = w'` |

★★`a ≠ 0` は `algebraMap a = (u : K)` が**単元**だから出る。
★★★したがって `tateXpair a w (a·w)`・`tateYpair a w (a·w)` は
**類 `c` だけの関数**である(証明項は `Prop` なので `rfl` で閉じる)。

## ★★残っているもの

(b) の厳密化(`I^n` 法の帰納、`adicSum_mul` が道具)、(c) 準同型性、
(d) の後半(核がちょうど `q^ℤ` であること = 単射性)、(e) 全射性。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `normRep` ほか | ★★★基本領域の代表元を関数にする |
| `eq_normRep` | ★★★★同じ類の正規化元は代表元に一致する |
| `exists_pair_of_class` | ★★★★どの類にも対 `(a, w)` が取れる |
| `pair_eq_of_same_class` | ★★★★★**対は類だけで決まる** |
| `tateXpair_eq_of_same_class` / `tateYpair_...` | ★★★★★★**座標は類だけで決まる** |
-/

namespace ABC3.Found.GaloisRep

open QuotientGroup

/-! ## ★★★基本領域の代表元を関数にする -/

section Rep

variable {K : Type} [Field K]

/-- ★★★基本領域の代表元を**関数として**取る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def normRep (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) : Kˣ :=
  (exists_unique_normalized_rep v q hq c).choose

theorem normRep_spec (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) :
    QuotientGroup.mk (normRep v q hq c) = c ∧ 0 ≤ vAdd v (normRep v q hq c)
      ∧ vAdd v (normRep v q hq c) < vAdd v q :=
  (exists_unique_normalized_rep v q hq c).choose_spec.1

theorem normRep_mk (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) : QuotientGroup.mk (normRep v q hq c) = c :=
  (normRep_spec v q hq c).1

theorem normRep_nonneg (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) : 0 ≤ vAdd v (normRep v q hq c) :=
  (normRep_spec v q hq c).2.1

theorem normRep_lt (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) : vAdd v (normRep v q hq c) < vAdd v q :=
  (normRep_spec v q hq c).2.2

/-- ★★★★**同じ類の正規化元は代表元に一致する**。 -/
theorem eq_normRep (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (c : Kˣ ⧸ Subgroup.zpowers q) (u : Kˣ) (hc : QuotientGroup.mk u = c)
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v q) :
    u = normRep v q hq c :=
  (exists_unique_normalized_rep v q hq c).choose_spec.2 u ⟨hc, h0, h1⟩

theorem normRep_eq_self (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q)
    (u : Kˣ) (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v q) :
    normRep v q hq (QuotientGroup.mk u) = u :=
  (eq_normRep v q hq _ u rfl h0 h1).symm

end Rep

/-! ## ★★★★★対は類だけで決まる -/

section Pair

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★**どの類にも対 `(a, w)` が取れる**。 -/
theorem exists_pair_of_class (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    (hmem : ∀ x : Kˣ, 0 < vAdd v x → ∃ y ∈ I, algebraMap R K y = (x : K))
    (hmem0 : ∀ x : Kˣ, 0 ≤ vAdd v x → ∃ y : R, algebraMap R K y = (x : K))
    (c : Kˣ ⧸ Subgroup.zpowers Q) :
    ∃ a w : R, algebraMap R K a = (normRep v Q hQ c : K) ∧ w ∈ I ∧
      algebraMap R K (a * w) = (Q : K) := by
  obtain ⟨a, w, ha, hwI, _, haw⟩ :=
    exists_pair_of_normalized (R := R) (I := I) v Q (normRep v Q hQ c)
      (normRep_nonneg v Q hQ c) (normRep_lt v Q hQ c) hmem hmem0
  exact ⟨a, w, ha, hwI, haw⟩

/-- ★★★**対は一意である**——`R → K` が単射で `R` が整域なら。 -/
theorem pair_unique [IsDomain R] (hinj : Function.Injective (algebraMap R K))
    {a w a' w' : R} (ha : algebraMap R K a = algebraMap R K a')
    (hq : algebraMap R K (a * w) = algebraMap R K (a' * w')) (ha0 : a ≠ 0) :
    a = a' ∧ w = w' := by
  have haa : a = a' := hinj ha
  subst haa
  have h : a * w = a * w' := hinj hq
  exact ⟨rfl, mul_left_cancel₀ ha0 h⟩

/-- ★★★★★**同じ類から来る 2 つの対は一致する**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`a ≠ 0` は `algebraMap a = (u : K)` が単元であることから出る。 -/
theorem pair_eq_of_same_class [IsDomain R] (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {a w a' w' : R} {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (hqa : algebraMap R K (a * w) = (Q : K)) (hqa' : algebraMap R K (a' * w') = (Q : K))
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    a = a' ∧ w = w' := by
  have huu : u = u' :=
    (eq_normRep v Q hQ _ u hcls h0 h1).trans (eq_normRep v Q hQ _ u' rfl h0' h1').symm
  have ha : algebraMap R K a = algebraMap R K a' := by rw [hu, hu', huu]
  have ha0 : a ≠ 0 := by
    intro h
    rw [h, map_zero] at hu
    exact (Units.ne_zero u) hu.symm
  exact pair_unique hinj ha (by rw [hqa, hqa']) ha0

/-- ★★★★★★**`X` は類だけで決まる**。 -/
theorem tateXpair_eq_of_same_class [IsDomain R] [IsAdicComplete I R]
    (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {a w a' w' : R} (hw : w ∈ I) (hw' : w' ∈ I) {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (hqa : algebraMap R K (a * w) = (Q : K)) (hqa' : algebraMap R K (a' * w') = (Q : K))
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    tateXpair (I := I) a w (a * w) (Ideal.mul_mem_left _ _ hw)
      = tateXpair (I := I) a' w' (a' * w') (Ideal.mul_mem_left _ _ hw') := by
  obtain ⟨rfl, rfl⟩ := pair_eq_of_same_class hinj v Q hQ hu hu' hqa hqa' h0 h1 h0' h1' hcls
  rfl

/-- ★★★★★★**`Y` は類だけで決まる**。 -/
theorem tateYpair_eq_of_same_class [IsDomain R] [IsAdicComplete I R]
    (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {a w a' w' : R} (hw : w ∈ I) (hw' : w' ∈ I) {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (hqa : algebraMap R K (a * w) = (Q : K)) (hqa' : algebraMap R K (a' * w') = (Q : K))
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    tateYpair (I := I) a w (a * w) (Ideal.mul_mem_left _ _ hw)
      = tateYpair (I := I) a' w' (a' * w') (Ideal.mul_mem_left _ _ hw') := by
  obtain ⟨rfl, rfl⟩ := pair_eq_of_same_class hinj v Q hQ hu hu' hqa hqa' h0 h1 h0' h1' hcls
  rfl

end Pair

/-! ## ★出典の紐付け(`.src`) -/

def normRep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——基本領域の代表元)",
    sectionId := "genell-def-3-3" }

def pair_eq_of_same_class.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——対 (a, w) が類だけで決まること)",
    sectionId := "genell-def-3-3" }

def tateXpair_eq_of_same_class.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——座標が K^x/q^Z 上で定まること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep

/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ClearDenominator
import ABC3.Meta.Claim

/-!
# ★★★★★★重なりの一致は冪で強制できる —— 段 E3a-3 の心臓（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 貼り合わせの障害を測る

`ClearDenominator.lean`（段 E3a の局所版）は「アフィン開 `U` の上では
`g ∈ Γ(X, X_s ⊓ U)` の分母を `f = trivValue s` の冪で払える」と言う。
★★しかし**それだけでは貼り合わない**——2 つのチャート `U_j, U_k` で得た `a_j, a_k` は
`U_j ⊓ U_k ⊓ X_s` の上でしか一致していないからである。

★★★**本ファイルがその障害を潰す**:

> `a` と `b` が `X.basicOpen f` の上で一致するなら、**ある `n` で `f^n·a = f^n·b`**。

★すなわち**もう一段 `s` の冪を掛ければ重なりでも一致する**。
これが [Stacks] Lemma 01PW の証明の第 2 段であり、原文が
「[some positive tensor power of]」という括弧に畳んでいる当のものである。

## ★★★★機構 —— 局所化の核だけ

    `IsLocalization.map_eq_zero_iff (Submonoid.powers f)`
      : `algebraMap R S a = 0 ↔ ∃ m ∈ powers f, m * a = 0`

★`a - b` に当てるだけである。`IsAffineOpen.isLocalization_basicOpen` が
`Γ(X, D(f))` を `Γ(X, U)` の `f` による局所化だと言うので、そこへ渡す。

## ★★★指数を揃える段（`exists_common_pow`）

有限被覆の重なりは有限個なので、各々の `n_{jk}` の**最大値**を取ればよい。
★そのために「単調な述語は有限個まとめて満たせる」という一般補題を置いた
——`Finset.sup` と `Finset.le_sup` だけである。

★★単調性は `pow_mul_eq_of_le`（`f^n·a = f^n·b` なら `f^m·a = f^m·b`、`n ≤ m`）が与える。

## ★残っている段（明示）

★★**貼り合わせそのもの**——`τ_j ∈ (M^{⊗N})(U_j)` が重なりで一致することを示したあと、
**層の条件で大域切断を取り出す**段は本ファイルに無い。
★§9-824 のとおり、それは `sheafify (M^{⊗N})` の側で行う（前層のままでは貼れない）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★橋渡し —— `algebraMap` は制限射そのものである -/

/-- ★`Γ(X, U) → Γ(X, D(f))` の `algebraMap` は前層の制限射そのものである（`rfl`）。

★★これで `ClearDenominator.lean` の主張（`algebraMap` で書いてある）を
貼り合わせの側（制限射で書く）へそのまま渡せる。 -/
theorem algebraMap_basicOpen_eq_res (U : X.Opens) (f : Γ(X, U)) (a : Γ(X, U)) :
    algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a
      = X.presheaf.map (homOfLE (X.basicOpen_le f)).op a := rfl

/-! ## ★★★★★★局所化の核 —— 消えるものは冪で消える -/

/-- ★★★★★★**`D(f)` の上で消える切断は `f` の冪を掛ければ消える**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★機構は `IsLocalization.map_eq_zero_iff (Submonoid.powers f)` だけである
——`IsAffineOpen.isLocalization_basicOpen` が局所化であることを与える。 -/
theorem exists_pow_mul_eq_zero (U : X.Opens) (hU : IsAffineOpen U) (f : Γ(X, U))
    (a : Γ(X, U))
    (h : algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a = 0) :
    ∃ n : ℕ, f ^ n * a = 0 := by
  haveI := hU.isLocalization_basicOpen f
  rw [IsLocalization.map_eq_zero_iff (Submonoid.powers f)] at h
  obtain ⟨⟨m, n, hn⟩, hm⟩ := h
  exact ⟨n, by simpa [hn] using hm⟩

/-- ★★★★★★★★**重なりの一致は冪で強制できる** —— 段 E3a-3 の心臓。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `a|_{D(f)} = b|_{D(f)}` ⟹ `∃ n, f^n·a = f^n·b`

★★これが貼り合わせの障害を潰す段である——2 つのチャートで分母を払って得た
`a_j, a_k` は `U_j ⊓ U_k ⊓ X_s` の上でしか一致しないが、
**もう一段 `s` の冪を掛ければ `U_j ⊓ U_k` 全体で一致する**。

★★★原文が「[some positive tensor power of]」という括弧に畳んでいるのは
まさにこの段である（分母を払う段と、重なりを合わせる段の 2 回、冪が要る）。 -/
theorem exists_pow_mul_eq_of_res_eq (U : X.Opens) (hU : IsAffineOpen U) (f : Γ(X, U))
    (a b : Γ(X, U))
    (h : algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a
       = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) b) :
    ∃ n : ℕ, f ^ n * a = f ^ n * b := by
  obtain ⟨n, hn⟩ := exists_pow_mul_eq_zero U hU f (a - b) (by rw [map_sub, h, sub_self])
  exact ⟨n, by rw [← sub_eq_zero, ← mul_sub]; exact hn⟩

/-! ## ★★★★★指数を揃える -/

/-- ★★**冪は上げてよい**（消滅の形）。 -/
theorem pow_mul_eq_zero_of_le (U : X.Opens) (f a : Γ(X, U)) {n : ℕ} (hn : f ^ n * a = 0)
    {m : ℕ} (hm : n ≤ m) : f ^ m * a = 0 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hm
  rw [pow_add, mul_assoc, mul_comm (f ^ k), ← mul_assoc, hn, zero_mul]

/-- ★★★**冪は上げてよい**（一致の形）——これが「最大値で揃える」を許す。 -/
theorem pow_mul_eq_of_le (U : X.Opens) (f a b : Γ(X, U)) {n : ℕ} (hn : f ^ n * a = f ^ n * b)
    {m : ℕ} (hm : n ≤ m) : f ^ m * a = f ^ m * b := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hm
  rw [pow_add, mul_assoc, mul_comm (f ^ k), ← mul_assoc, hn]
  ring

open scoped Classical in
/-- ★★★★★**単調な述語は有限個まとめて満たせる** —— 「指数を最大値で揃える」段。

★有限被覆の重なりは有限個なので、各々の `n_{jk}` の最大値を取ればよい。
★★機構は `Finset.sup` と `Finset.le_sup` だけである。
単調性は `pow_mul_eq_of_le` / `pow_mul_eq_zero_of_le` が与える。 -/
theorem exists_common_pow {ι : Type} (I : Finset ι) (P : ι → ℕ → Prop)
    (mono : ∀ i n m, n ≤ m → P i n → P i m) (h : ∀ i ∈ I, ∃ n, P i n) :
    ∃ N : ℕ, ∀ i ∈ I, P i N := by
  choose! n hn using h
  exact ⟨I.sup n, fun i hi => mono i (n i) _ (Finset.le_sup hi) (hn i hi)⟩

/-- ★★★★★★★**有限個の重なりで一斉に一致させる** —— 段 E3a-3 の第 2 段そのもの。

各 `i ∈ I` について、アフィン開 `U i` の上の `a i`・`b i` が `D(f i)` で一致するなら、
★**すべての `i` で同時に効く単一の指数 `N`** がある。 -/
theorem exists_common_pow_mul_eq {ι : Type} (I : Finset ι) (U : ι → X.Opens)
    (hU : ∀ i ∈ I, IsAffineOpen (U i)) (f a b : ∀ i, Γ(X, U i))
    (h : ∀ i ∈ I, algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) (a i)
       = algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) (b i)) :
    ∃ N : ℕ, ∀ i ∈ I, (f i) ^ N * a i = (f i) ^ N * b i := by
  refine exists_common_pow I (fun i n => (f i) ^ n * a i = (f i) ^ n * b i)
    (fun i n m hnm hn => pow_mul_eq_of_le (U i) (f i) (a i) (b i) hn hnm) ?_
  exact fun i hi => exists_pow_mul_eq_of_res_eq (U i) (hU i hi) (f i) (a i) (b i) (h i hi)

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_mul_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(D(f) の上で消える切断は f の冪で消える)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_of_res_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの一致は f の冪で強制できる)",
    sectionId := "genell-prop-1-4" }

def exists_common_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(指数を最大値で揃える)",
    sectionId := "genell-prop-1-4" }

def exists_common_pow_mul_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(有限個の重なりで一斉に一致させる)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_of_res_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsLocalization.map_eq_zero_iff(局所化の核)"
      (.inMathlib "IsLocalization.map_eq_zero_iff") 7,
    .citation "[mathlib]" "IsAffineOpen.isLocalization_basicOpen"
      (.inMathlib "AlgebraicGeometry.IsAffineOpen.isLocalization_basicOpen") 7,
    .citation "[Stacks]" "Lemma 01PW の第 2 段(重なりを合わせる段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★原文の「[some positive tensor power of]」という括弧には**冪が 2 回**畳まれている" ++
       "——(1) 分母を払う段(ClearDenominator.lean)と (2) 重なりを合わせる段(本ファイル)。" ++
       "★★どちらも局所化の性質だけで出るが、**別々の補題である**") 7,
    .implicitStep
      ("★★**貼り合わせそのものは本ファイルに無い**。重なりで一致させたあと" ++
       "層の条件で大域切断を取り出す段が残る。§9-824 のとおりそれは " ++
       "sheafify (M^{⊗N}) の側で行う——前層のままでは貼れない") 8 ]

end ABC3.Found.GenEll

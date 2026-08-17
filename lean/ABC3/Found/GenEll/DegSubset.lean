import ABC3.Found.GenEll.UPoint

/-!
# [GenEll] Example 1.3, (i) —— **次数で切った点の部分集合と Galois 有限性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> (i) If d ∈ N {∞}, then we shall denote by

## ★★★これは定理ではなく**定義**である

原文の (i) は 4 つの記号を導入するだけである:

| 記号 | 意味 |
|---|---|
| `X(ℚ̄)^{≤d}` | `[F:ℚ] ≤ d` なる数体 `F` について `X(F)` を合併したもの |
| `X(ℚ̄)^{=d}` | `X(ℚ̄)^{≤d} \ X(ℚ̄)^{≤d-1}` |
| `E^{≤d}`, `E^{=d}` | 部分集合 `E` との交わり |
| Galois-finite | 各 `E^{≤d}`(`d` は正整数)が有限 |

★★`d` は **`ℕ ∪ {∞}`** を走る。本ファイルは `ℕ∞` で受ける
(原文の `X(ℚ̄)^{≤∞} = X(ℚ̄)` が `leDeg_top` である)。

## ★★★原文が指している「最小定義体」との接続

原文は (i) の直後に
> [cf. the discussion in Definition 1.5, (i), below of "minimal fields of definition"].

と書いている。★これは `X(ℚ̄)^{=d}` の読み方の指示である——
`x ∈ X(ℚ̄)^{=d}` とは「`x` の定義体の次数の**最小値が `d`**」ということである。
★★それを `mem_eqDeg_iff` として証明した。**差集合の定義から出る**。

## ★点の型は `UPoint`

`UPoint X D` は底変換で同一視した代数的点の型(`UPoint.lean`)。
★★★**「次数」は代表元ごとに違う**——だから `≤ d` は
**「そういう代表元が存在する」**として定義するのが正しい。
★それが `Quot` の上で well-defined であることは、存在量化なので自動である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable {X : Scheme.{0}} {D : X.IdealSheafData}

/-! ## ★定義体の次数 -/

/-- ★**代表元の定義体の次数** `[F : ℚ]`。 -/
noncomputable def AlgPointOff.deg (p : AlgPointOff X D) : ℕ :=
  letI := p.instField
  letI := p.instNF
  Module.finrank ℚ p.fld

/-- ★数体の次数は正である。★★これが `leDeg_zero` を出す。 -/
theorem AlgPointOff.deg_pos (p : AlgPointOff X D) : 0 < p.deg := by
  letI := p.instField
  letI := p.instNF
  exact Module.finrank_pos

@[simp] theorem deg_algPointOff (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X) (h : pullbackIdeal F D xF ≠ 0) :
    (algPointOff F xF h).deg = Module.finrank ℚ F := rfl

/-! ## ★★★`X(ℚ̄)^{≤d}` -/

/-- ★★★**[GenEll] Example 1.3, (i)** —— `X(ℚ̄)^{≤d}`。

原文 (GenEll p.5):
> (i) If d ∈ N {∞}, then we shall denote by

    `X(ℚ̄)^{≤d}` = `[F:ℚ] ≤ d` なる `F` についての `X(F)` の合併

★★★**「合併」を「そういう代表元が存在する」と読む**のが、
商 `UPoint` の上での正しい形である。 -/
def leDeg (X : Scheme.{0}) (D : X.IdealSheafData) (d : ℕ∞) : Set (UPoint X D) :=
  {x | ∃ p : AlgPointOff X D, UPoint.mk p = x ∧ (p.deg : ℕ∞) ≤ d}

theorem mem_leDeg_mk (p : AlgPointOff X D) {d : ℕ∞} (h : (p.deg : ℕ∞) ≤ d) :
    UPoint.mk p ∈ leDeg X D d := ⟨p, rfl, h⟩

theorem mem_leDeg_self (p : AlgPointOff X D) : UPoint.mk p ∈ leDeg X D p.deg :=
  mem_leDeg_mk p le_rfl

/-- ★`d` について単調。 -/
theorem leDeg_mono {d e : ℕ∞} (h : d ≤ e) : leDeg X D d ⊆ leDeg X D e :=
  fun _ hx => hx.imp fun _ hp => ⟨hp.1, hp.2.trans h⟩

/-- ★★★**`X(ℚ̄)^{≤∞} = X(ℚ̄)`**(原文が明記している)。 -/
@[simp] theorem leDeg_top : leDeg X D ⊤ = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  obtain ⟨p, hp⟩ := Quot.exists_rep x
  exact ⟨p, hp, le_top⟩

/-- ★**`X(ℚ̄)^{≤0} = ∅`**——数体の次数は正だからである。

★★これがあるので Galois 有限性の「`d` は正整数を走る」という但し書きは
**制限になっていない**(`galoisFinite_iff`)。 -/
@[simp] theorem leDeg_zero : leDeg X D 0 = ∅ := by
  ext x
  simp only [Set.mem_empty_iff_false, iff_false]
  rintro ⟨p, -, hd⟩
  exact absurd (le_antisymm (by exact_mod_cast hd) (Nat.zero_le _)) p.deg_pos.ne'

/-- ★代表元があれば、その次数以上のところには入っている。 -/
theorem notMem_leDeg_iff (x : UPoint X D) (d : ℕ∞) :
    x ∉ leDeg X D d ↔ ∀ p : AlgPointOff X D, UPoint.mk p = x → d < (p.deg : ℕ∞) := by
  simp only [leDeg, Set.mem_setOf_eq, not_exists, not_and, not_le]

/-! ## ★★`X(ℚ̄)^{=d}` -/

/-- ★★**[GenEll] Example 1.3, (i)** —— `X(ℚ̄)^{=d} = X(ℚ̄)^{≤d} \ X(ℚ̄)^{≤d-1}`。

原文 (GenEll p.5):
> (i) If d ∈ N {∞}, then we shall denote by
-/
def eqDeg (X : Scheme.{0}) (D : X.IdealSheafData) (d : ℕ∞) : Set (UPoint X D) :=
  leDeg X D d \ leDeg X D (d - 1)

theorem eqDeg_subset_leDeg (d : ℕ∞) : eqDeg X D d ⊆ leDeg X D d :=
  Set.sdiff_subset

/-- ★★`X(ℚ̄)^{≤d}` は 1 段下と `X(ℚ̄)^{=d}` に分かれる。 -/
theorem leDeg_eq_union (d : ℕ∞) :
    leDeg X D d = leDeg X D (d - 1) ∪ eqDeg X D d := by
  rw [eqDeg, Set.union_sdiff_self]
  exact (Set.union_eq_self_of_subset_left (leDeg_mono tsub_le_self)).symm

/-- ★`X(ℚ̄)^{=∞} = ∅`——`∞ - 1 = ∞` だからである。

★★原文が `d ∈ ℕ ∪ {∞}` と書いているとおりに読むと、こうなる。 -/
@[simp] theorem eqDeg_top : eqDeg X D ⊤ = ∅ := by
  simp [eqDeg]

/-- ★★★**原文が `Definition 1.5, (i)` を参照している理由**。

原文 (GenEll p.5):
> [cf. the discussion in Definition 1.5, (i), below of "minimal fields of definition"].

★★`x ∈ X(ℚ̄)^{=d}` とは、**`x` の定義体の次数の最小値が `d`** ということである
——すなわち「次数 `d` の代表元があり、かつどの代表元も次数 `d` 以上」。

★★★**差集合の定義から出る。**原文の但し書きはこの読み方の指示である。 -/
theorem mem_eqDeg_iff (x : UPoint X D) (d : ℕ) :
    x ∈ eqDeg X D (d + 1) ↔
      (∃ p : AlgPointOff X D, UPoint.mk p = x ∧ p.deg = d + 1) ∧
      (∀ p : AlgPointOff X D, UPoint.mk p = x → d + 1 ≤ p.deg) := by
  have hsub : ((d : ℕ∞) + 1) - 1 = (d : ℕ∞) := by
    rw [← Nat.cast_one (R := ℕ∞), ← Nat.cast_add, ← ENat.coe_sub]
    simp
  constructor
  · rintro ⟨⟨p, hp, hle⟩, hout⟩
    rw [hsub] at hout
    have key : ∀ q : AlgPointOff X D, UPoint.mk q = x → d + 1 ≤ q.deg := by
      intro q hq
      have := (notMem_leDeg_iff x (d : ℕ∞)).1 hout q hq
      exact_mod_cast Order.add_one_le_of_lt (by exact_mod_cast this)
    refine ⟨⟨p, hp, le_antisymm ?_ (key p hp)⟩, key⟩
    exact_mod_cast hle
  · rintro ⟨⟨p, hp, hd⟩, hall⟩
    refine ⟨⟨p, hp, by rw [hd]; exact_mod_cast le_rfl⟩, ?_⟩
    rw [hsub]
    refine (notMem_leDeg_iff x (d : ℕ∞)).2 fun q hq => ?_
    have := hall q hq
    exact_mod_cast Nat.lt_of_succ_le this

/-! ## ★★部分集合の `E^{≤d}` / `E^{=d}` -/

/-- ★**[GenEll] Example 1.3, (i)** —— `E^{≤d} = E ∩ X(ℚ̄)^{≤d}`。 -/
def subLeDeg (E : Set (UPoint X D)) (d : ℕ∞) : Set (UPoint X D) := E ∩ leDeg X D d

/-- ★**[GenEll] Example 1.3, (i)** —— `E^{=d} = E ∩ X(ℚ̄)^{=d}`。 -/
def subEqDeg (E : Set (UPoint X D)) (d : ℕ∞) : Set (UPoint X D) := E ∩ eqDeg X D d

theorem subLeDeg_mono_set {E E' : Set (UPoint X D)} (h : E' ⊆ E) (d : ℕ∞) :
    subLeDeg E' d ⊆ subLeDeg E d :=
  Set.inter_subset_inter_left _ h

theorem subEqDeg_subset_subLeDeg (E : Set (UPoint X D)) (d : ℕ∞) :
    subEqDeg E d ⊆ subLeDeg E d :=
  Set.inter_subset_inter_right _ (eqDeg_subset_leDeg d)

@[simp] theorem subLeDeg_top (E : Set (UPoint X D)) : subLeDeg E ⊤ = E := by
  simp [subLeDeg]

/-! ## ★★★Galois 有限性 -/

/-- ★★★**[GenEll] Example 1.3, (i)** —— `E` が **Galois-finite** であること。

原文 (GenEll p.5):
> over the positive integers] is finite, then we shall say that E is Galois-finite.

★原文どおり **`d` は正整数**を走る。★★`X(ℚ̄)^{≤0} = ∅` なので
`d = 0` を許しても同じである(`galoisFinite_iff`)。 -/
def GaloisFinite (E : Set (UPoint X D)) : Prop :=
  ∀ d : ℕ, 0 < d → (subLeDeg E (d : ℕ∞)).Finite

/-- ★**「正整数」という但し書きは制限になっていない。** -/
theorem galoisFinite_iff (E : Set (UPoint X D)) :
    GaloisFinite E ↔ ∀ d : ℕ, (subLeDeg E (d : ℕ∞)).Finite := by
  refine ⟨fun h d => ?_, fun h d _ => h d⟩
  rcases Nat.eq_zero_or_pos d with rfl | hd
  · simp [subLeDeg]
  · exact h d hd

theorem GaloisFinite.mono {E E' : Set (UPoint X D)} (h : E' ⊆ E) (hE : GaloisFinite E) :
    GaloisFinite E' :=
  fun d hd => (hE d hd).subset (subLeDeg_mono_set h _)

theorem galoisFinite_of_finite {E : Set (UPoint X D)} (h : E.Finite) : GaloisFinite E :=
  fun _ _ => h.subset Set.inter_subset_left

@[simp] theorem galoisFinite_empty : GaloisFinite (∅ : Set (UPoint X D)) :=
  galoisFinite_of_finite Set.finite_empty

theorem GaloisFinite.union {E E' : Set (UPoint X D)}
    (hE : GaloisFinite E) (hE' : GaloisFinite E') : GaloisFinite (E ∪ E') := by
  intro d hd
  have : subLeDeg (E ∪ E') (d : ℕ∞) = subLeDeg E d ∪ subLeDeg E' d :=
    Set.union_inter_distrib_right _ _ _
  rw [this]
  exact (hE d hd).union (hE' d hd)

theorem GaloisFinite.inter_left {E : Set (UPoint X D)} (hE : GaloisFinite E)
    (E' : Set (UPoint X D)) : GaloisFinite (E ∩ E') :=
  hE.mono Set.inter_subset_left

/-- ★★Galois 有限なら各 `E^{=d}` も有限である。 -/
theorem GaloisFinite.subEqDeg_finite {E : Set (UPoint X D)} (hE : GaloisFinite E)
    (d : ℕ) (hd : 0 < d) : (subEqDeg E (d : ℕ∞)).Finite :=
  (hE d hd).subset (subEqDeg_subset_subLeDeg E _)

/-! ## ★出典の紐付け(`.src`)

★★★**条なしである。**原文の `Example 1.3, (i)` は
4 つの記号(`X(ℚ̄)^{≤d}` / `X(ℚ̄)^{=d}` / `E^{≤d}`,`E^{=d}` / Galois-finite)を
導入するだけであり、**そのすべてを定義した**。

★原文が付けている但し書き 2 つ——`X(ℚ̄)^{≤∞} = X(ℚ̄)` と
「最小定義体」への参照——も定理として取った
(`leDeg_top` / `mem_eqDeg_iff`)。 -/

def leDeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

def GaloisFinite.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

end ABC3.Found.GenEll

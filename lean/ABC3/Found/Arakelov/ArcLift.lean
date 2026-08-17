import ABC3.Found.Arakelov.ArcLandsIn
import ABC3.Found.Arakelov.ArcOpenImmersion

/-!
# Arakelov (C1) 最後の段 —— **開埋め込みへの持ち上げ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★C1 の残り 1 向きの最後の段

`topology_openImmersion` の残り 1 向きは「`(· ≫ f)` が開写像」である。
★手順(`genell-goal.md` §9-23)の 3 段のうち、

| 段 | 状態 |
|---|---|
| 1. chart 逆像を評価に帰着 | ✅ `mem_basicOpen_base_iff` |
| 2. 「点が開集合に落ちる」が開 | ✅ `isOpen_landsIn` |
| 3. **逆写像 `ψ` の構成と連続性** | ★★★**本ファイル** |

## ★★★像への所属は「点が像に落ちる」ことと同値である

`q : Spec ℂ ⟶ Y` が `(· ≫ f)` の像に入る ⟺ `q` の像の点が `f` の像に入る。

★★★**これが「像が開である」の内実**である——右辺は
`isOpen_landsIn` により**開条件**だからである。

★機構は mathlib の `IsOpenImmersion.lift`
(像が入っていれば持ち上がる)と `lift_fac` / `lift_uniq`。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★★像への所属の特徴づけ -/

/-- ★★★**`(· ≫ f)` の像に入る ⟺ 点が `f` の像に落ちる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが「像が開である」の内実である**——右辺は `isOpen_landsIn` により開条件。

★機構は `IsOpenImmersion.lift`(⟸)と合成の像の計算(⟹)。 -/
theorem mem_range_comp_iff {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (q : Spec (CommRingCat.of ℂ) ⟶ Y) :
    q ∈ Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f)
      ↔ Set.range q.base ⊆ Set.range f.base := by
  constructor
  · rintro ⟨p, rfl⟩
    intro y hy
    obtain ⟨x, rfl⟩ := hy
    exact ⟨p.base x, rfl⟩
  · intro h
    exact ⟨IsOpenImmersion.lift f q h, IsOpenImmersion.lift_fac f q h⟩

/-- ★★**像に入る条件は `Spec ℂ` の唯一の点だけで判定できる**。 -/
theorem mem_range_comp_iff_base_default {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (q : Spec (CommRingCat.of ℂ) ⟶ Y) :
    q ∈ Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f)
      ↔ q.base default ∈ Set.range f.base := by
  rw [mem_range_comp_iff]
  constructor
  · intro h; exact h ⟨default, rfl⟩
  · intro h y hy
    obtain ⟨x, rfl⟩ := hy
    have : x = default := Subsingleton.elim _ _
    rw [this]
    exact h

/-! ## ★★持ち上げの一意性 -/

/-- ★**持ち上げは一意である**(開埋め込みはモノだから)。 -/
theorem lift_eq_of_comp_eq {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    {q : Spec (CommRingCat.of ℂ) ⟶ Y} (h : Set.range q.base ⊆ Set.range f.base)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ≫ f = q) :
    p = IsOpenImmersion.lift f q h :=
  IsOpenImmersion.lift_uniq f q h p hp

/-! ## ★★★持ち上げ写像(点の条件から) -/

/-- ★**唯一の点での条件から、像全体の包含が出る**(`Spec ℂ` は一点だから)。 -/
theorem range_base_subset_of_mem {X Y : Scheme.{0}} (f : X ⟶ Y)
    (q : Spec (CommRingCat.of ℂ) ⟶ Y) (h : q.base default ∈ Set.range f.base) :
    Set.range q.base ⊆ Set.range f.base := by
  intro y hy
  obtain ⟨x, rfl⟩ := hy
  have hx : x = default := Subsingleton.elim _ _
  rw [hx]
  exact h

/-- ★★★**開埋め込みへの持ち上げ**——像の点が `f` の像に入っていれば持ち上がる。 -/
noncomputable def openLift {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (q : Spec (CommRingCat.of ℂ) ⟶ Y) (h : q.base default ∈ Set.range f.base) :
    Spec (CommRingCat.of ℂ) ⟶ X :=
  IsOpenImmersion.lift f q (range_base_subset_of_mem f q h)

/-- ★★**持ち上げてから合成するともとに戻る**。 -/
@[simp] theorem openLift_comp {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (q : Spec (CommRingCat.of ℂ) ⟶ Y) (h : q.base default ∈ Set.range f.base) :
    openLift f q h ≫ f = q :=
  IsOpenImmersion.lift_fac f q _

/-- ★★**逆に、合成してから持ち上げるともとに戻る**(単射性から)。 -/
theorem openLift_comp_self {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : (p ≫ f).base default ∈ Set.range f.base) :
    openLift f (p ≫ f) h = p := by
  exact comp_openImmersion_injective f (openLift_comp f (p ≫ f) h)

/-! ## ★出典の紐付け(`.src`) -/

def mem_range_comp_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——像への所属が点の条件で決まること)",
    sectionId := "genell-def-1-1-i" }

def openLift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みへの持ち上げ写像)",
    sectionId := "genell-def-1-1-i" }

def mem_range_comp_iff_base_default.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——像への所属が唯一の点で判定できること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov

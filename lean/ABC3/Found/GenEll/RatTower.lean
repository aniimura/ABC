/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Localization.FractionRing
import Mathlib.Data.Nat.Factorial.BigOperators
import Mathlib.Algebra.Category.Ring.Colimits
import Mathlib.CategoryTheory.Functor.OfSequence
import ABC3.Meta.Claim

/-!
# [GenEll] Remark 1.5.1 —— **`ℤ[Σ⁻¹]` の増大列**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★これは spreading out の第一石である

`Remark 1.5.1` の後半（`log-cond_D` の BD-類が `(X_ℚ, D_ℚ)` にしか依らないこと）の
証明は、**生成ファイバーの同型が有限個の素数 Σ を除いて延びる**ことを使う。
★その台となるのが `ℤ[Σ⁻¹]` の増大列と、その合併が `ℚ` であることである。

## ★★★★添字を `Σ` でなく `n!` にした

素数の有限集合 `Σ` で添字を取ると `ℤ[Σ⁻¹] ⊆ ℤ[Σ'⁻¹]` を言うのに `Σ ⊆ Σ'` が要り、
**整除性を毎回作る**ことになる。★そこで `ℤ[1/n!]` で添字を取る:

    n ≤ m  ⟹  n! ∣ m!  （`Nat.factorial_dvd_factorial`）

——★★**順序がそのまま整除性になる**ので遷移が自動で入る。
★★★どの有理数 `q` も `q.den ∣ q.den!` なので `ℤ[1/q.den!]` に入る（`mem_ratTower_den`）。

## ★★★★★到達点

| 主張 | 宣言 |
|---|---|
| 増大列 | `ratTower_mono` |
| 汲み尽くし（各有理数） | `mem_ratTower_den` |
| 有限集合はある段に収まる | ★`exists_ratTower_of_finset` |
| ★**アフィンの spreading out** | ★★`exists_ratTower_range_le` |

★★`exists_ratTower_range_le` —— **有限生成 ℤ-代数から `ℚ` への射は `ℤ[1/n!]` を経由する**。
これがアフィンな場合の spreading out そのものである。

## ★★★★★★mathlib 側の状況（2026-08-27 の実測）

★**`Mathlib/AlgebraicGeometry/AffineTransitionLimit.lean` が EGA IV §8 である**
（`ResearchPaper/mathlib-gap.json` の `spreading-out-over-base` を同日訂正）。

| 段 | 状態 |
|---|---|
| `lim D_i ⟶ X` が有限段を経由する | ★mathlib `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` |
| 単射側（Stacks 01ZC） | ★mathlib `Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType` |
| `Spec ℚ = lim_n Spec ℤ[1/n!]` の図式 | ★**本ファイルが台**（`IsLimit` の組み立ては未） |
| 同型であることの降下 | ★★上の 2 つから出る（両向きの射を取り、合成が恒等と生成ファイバーで一致 → 単射側で有限段でも一致） |

★★★**「同型の降下」に mathlib の欠落は無い**——射の降下と単射側から組める。
（当初は Stacks 01ZM が要ると見ていたが、両向きを取れば要らない。）
-/

namespace ABC3.Found.GenEll

open Nat

/-! ## ★★★★`ℤ[1/n!] ⊆ ℚ` -/

/-- ★**`ℤ[1/n!]`** —— `ℚ` の部分環。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★添字を素数の有限集合 `Σ` ではなく `n!` に取ることで、順序がそのまま整除性になる。 -/
def ratTower (n : ℕ) : Subring ℚ := Subring.closure {(((n ! : ℕ) : ℚ))⁻¹}

theorem intCast_mem_ratTower (n : ℕ) (a : ℤ) : (a : ℚ) ∈ ratTower n := intCast_mem _ a

theorem inv_factorial_mem_ratTower (n : ℕ) : (((n ! : ℕ) : ℚ))⁻¹ ∈ ratTower n :=
  Subring.subset_closure rfl

/-- ★★**増大列** —— `n ≤ m` なら `ℤ[1/n!] ⊆ ℤ[1/m!]`。

★`n! ∣ m!` なので `1/n! = (m!/n!)·(1/m!)` である。 -/
theorem ratTower_mono {n m : ℕ} (h : n ≤ m) : ratTower n ≤ ratTower m := by
  rw [ratTower, Subring.closure_le]
  rintro x rfl
  obtain ⟨c, hc⟩ : (n ! : ℕ) ∣ (m ! : ℕ) := Nat.factorial_dvd_factorial h
  have hcne : (c : ℚ) ≠ 0 := by
    intro h0
    have hc0 : c = 0 := by exact_mod_cast h0
    rw [hc0, mul_zero] at hc
    exact (Nat.factorial_ne_zero m) hc
  have hx : ((n ! : ℕ) : ℚ)⁻¹ = (c : ℚ) * ((m ! : ℕ) : ℚ)⁻¹ := by
    rw [hc]
    push_cast
    rw [mul_inv, ← mul_assoc, mul_comm (c:ℚ), mul_assoc, mul_inv_cancel₀ hcne, mul_one]
  simp only [SetLike.mem_coe]
  rw [hx]
  exact mul_mem (natCast_mem _ c) (inv_factorial_mem_ratTower m)

/-! ## ★★★★★★汲み尽くし -/

/-- ★★★**どの有理数もある段に入る** —— `q ∈ ℤ[1/q.den!]`。

★`q.den ∣ q.den!` なので `q = q.num·(q.den!/q.den)·(1/q.den!)` である。 -/
theorem mem_ratTower_den (q : ℚ) : q ∈ ratTower q.den := by
  obtain ⟨c, hc⟩ : (q.den : ℕ) ∣ (q.den ! : ℕ) := Nat.dvd_factorial q.pos le_rfl
  have hfne : ((q.den ! : ℕ) : ℚ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero q.den
  have key : q * ((q.den ! : ℕ) : ℚ) = (q.num : ℚ) * (c : ℚ) := by
    rw [hc]
    push_cast
    rw [← mul_assoc, Rat.mul_den_eq_num]
  have hq : (q.num : ℚ) * (c : ℚ) * ((q.den ! : ℕ) : ℚ)⁻¹ = q := by
    rw [← key, mul_assoc, mul_inv_cancel₀ hfne, mul_one]
  have hmem : (q.num : ℚ) * (c : ℚ) * ((q.den ! : ℕ) : ℚ)⁻¹ ∈ ratTower q.den :=
    mul_mem (mul_mem (intCast_mem_ratTower _ _) (natCast_mem _ c))
      (inv_factorial_mem_ratTower _)
  rwa [hq] at hmem

/-- ★★★★**有限集合はある段に収まる** —— 分母の最大値を取るだけ。

★★これが「有限個の素数 `Σ` を除いて」の中身である。 -/
theorem exists_ratTower_of_finset (s : Finset ℚ) :
    ∃ n : ℕ, ∀ q ∈ s, q ∈ ratTower n := by
  classical
  refine ⟨(s.image Rat.den).sup id, fun q hq => ?_⟩
  refine ratTower_mono ?_ (mem_ratTower_den q)
  exact Finset.le_sup (f := id) (Finset.mem_image_of_mem Rat.den hq)

/-! ## ★★★★★★★★アフィンの spreading out -/

/-- ★★★★★★★★**アフィンの spreading out** ——
有限生成 ℤ-代数から `ℚ` への射は `ℤ[1/n!]` を経由する。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★★**これがアフィンな場合の spreading out そのものである**——
生成元の像の分母を集めて共通の段を取るだけ。
★一般（非アフィン）の場合は mathlib の
`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`（EGA IV §8、
`AffineTransitionLimit.lean`）が受け持つが、そこへ渡す図式
`Spec ℚ = lim_n Spec ℤ[1/n!]` の台が本ファイルである。 -/
theorem exists_ratTower_range_le {S : Type} [CommRing S] (φ : S →+* ℚ)
    (s : Finset S) (hs : Subring.closure (s : Set S) = ⊤) :
    ∃ n : ℕ, ∀ x : S, φ x ∈ ratTower n := by
  classical
  obtain ⟨n, hn⟩ := exists_ratTower_of_finset (s.image φ)
  refine ⟨n, fun x => ?_⟩
  have hx : x ∈ Subring.closure (s : Set S) := by rw [hs]; trivial
  refine Subring.closure_induction (p := fun y _ => φ y ∈ ratTower n) ?_ ?_ ?_ ?_ ?_ ?_ hx
  · intro y hy
    exact hn (φ y) (Finset.mem_image_of_mem φ hy)
  · simp
  · simp
  · intro a b _ _ ha hb; rw [map_add]; exact add_mem ha hb
  · intro a _ ha; rw [map_neg]; exact neg_mem ha
  · intro a b _ _ ha hb; rw [map_mul]; exact mul_mem ha hb

/-! ## ★★★★★★★★余極限へ —— 図式と余錐

★★`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`（mathlib、EGA IV §8）へ
渡すには `IsLimit` の錐が要る。その台となる**環の側の図式と余錐**をここに置く。 -/

/-- ★★**有向性** —— どの 2 段も上に共通の段を持つ。 -/
theorem ratTower_directed : Directed (· ≤ ·) ratTower := by
  intro n m
  exact ⟨max n m, ratTower_mono (le_max_left n m), ratTower_mono (le_max_right n m)⟩

/-- ★★★★**合併は `ℚ` 全体** —— 汲み尽くしの束論的な形。 -/
theorem iSup_ratTower_eq_top : (⨆ n : ℕ, ratTower n) = ⊤ := by
  refine eq_top_iff.2 (fun q _ => ?_)
  exact le_iSup ratTower q.den (mem_ratTower_den q)

open CategoryTheory Limits

/-- ★★★**図式** `n ↦ ℤ[1/n!]`。

★`Functor.ofSequence` を使う——`n ⟶ n+1` の射だけ与えれば `map_id` / `map_comp` は
mathlib が持つ。★★**素朴に `Functor` の構造体を書くと核の判定が発散する**
（2026-08-27 に実測。`Subring.inclusion` の defeq が重い）。 -/
def ratTowerDiagram : ℕ ⥤ CommRingCat.{0} :=
  Functor.ofSequence (X := fun n => CommRingCat.of (ratTower n))
    (fun n => CommRingCat.ofHom (Subring.inclusion (ratTower_mono (Nat.le_succ n))))

/-- ★★★**余錐** —— 頂点は `ℚ`、成分は包含。

★`NatTrans.ofSequence` で自然性も `n ⟶ n+1` だけで済む。 -/
def ratTowerCocone : Cocone ratTowerDiagram where
  pt := CommRingCat.of ℚ
  ι := NatTrans.ofSequence (F := ratTowerDiagram) (G := (Functor.const ℕ).obj (CommRingCat.of ℚ))
    (fun n => CommRingCat.ofHom (ratTower n).subtype)
    (fun n => by
      simp only [ratTowerDiagram, Functor.ofSequence_map_homOfLE_succ,
        Functor.const_obj_map]
      apply CommRingCat.hom_ext; apply RingHom.ext; intro x
      rfl)

/-! ### ★★★★★★次の一歩 —— `IsColimit`

★この余錐が余極限であること（＝`ℚ = colim ℤ[1/n!]`）が次の段である。

★★**道は測ってある**（2026-08-27）: mathlib の
`Subalgebra.iSupLift`（`Mathlib/Algebra/Algebra/Subalgebra/Directed.lean`）が
「有向な部分代数の上限からの射を各段の射から作る」を与える。
`ratTower_directed` と `iSup_ratTower_eq_top` がその 2 つの入力である。

★★★ただし `iSupLift` は **`Subalgebra R A` 用で `Subring` 版が無い**。
`ratTower n` を `Subalgebra ℤ ℚ` として取り直すか、`toSubring` で橋を渡す必要がある。

★★★★そのあと `Spec`（`ΓSpec.adjunction` の右随伴）で `IsLimit` に移し、
`X ×_ℤ (−)`（`Over.mapPullbackAdj` の右随伴）で `X_ℚ = lim (X ×_ℤ ℤ[1/n!])` にする。 -/

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` の後半は
(1) 生成ファイバーの同型が `ℤ[Σ⁻¹]` へ延びる、(2) 因子 `D` も一緒に延びる、
(3) `Σ` の外で conductor が一致する、の 3 段からなり、
本ファイルは **(1) の台**（増大列と汲み尽くし、およびアフィンな場合）だけである。 -/

def ratTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(ℤ[Σ⁻¹] の増大列——添字を n! に取った)",
    sectionId := "genell-rem-1-5-1" }

def mem_ratTower_den.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(汲み尽くし——どの有理数もある段に入る)",
    sectionId := "genell-rem-1-5-1" }

def exists_ratTower_of_finset.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(有限集合はある段に収まる)",
    sectionId := "genell-rem-1-5-1" }

def exists_ratTower_range_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(アフィンの spreading out——非アフィンの図式は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_ratTower_range_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_ratTower_of_finset(有限集合はある段に収まる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_ratTower_of_finset") 9,
    .citation "[mathlib]" "Nat.factorial_dvd_factorial(n ≤ m ⟹ n! ∣ m! ——遷移射の出どころ)"
      (.inMathlib "Nat.factorial_dvd_factorial") 9,
    .implicitStep
      ("★★★Remark 1.5.1 の後半は (1) 生成ファイバーの同型が ℤ[Σ⁻¹] へ延びる、" ++
       "(2) 因子 D も一緒に延びる、(3) Σ の外で conductor が一致する、の 3 段。" ++
       "本ファイルは (1) の台(増大列と汲み尽くし、およびアフィンな場合)だけである") 9,
    .implicitStep
      ("★★2026-08-27 の実測: 非アフィンの段は mathlib の " ++
       "Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation" ++
       "(AffineTransitionLimit.lean、EGA IV §8)が受け持つ。そこへ渡す図式 " ++
       "Spec ℚ = lim_n Spec ℤ[1/n!] の IsLimit の組み立てが残っている") 9,
    .implicitStep
      ("★★★「同型であることの降下」に mathlib の欠落は無い——両向きの射を降ろし、" ++
       "合成が恒等と生成ファイバーで一致することを単射側" ++
       "(exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType)に流せば有限段で一致する。" ++
       "当初 Stacks 01ZM が要ると見ていたが要らない") 9,
    .implicitStep
      ("★添字を素数の有限集合 Σ ではなく n! に取った。順序 n ≤ m がそのまま整除性 " ++
       "n! ∣ m! になるので遷移射が自動で入る。原文の Σ との対応は、" ++
       "n! の素因数の集合を Σ と読むことで付く") 9 ]

end ABC3.Found.GenEll

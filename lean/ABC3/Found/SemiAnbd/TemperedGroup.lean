import Mathlib.Topology.Algebra.ContinuousMonoidHom
import Mathlib.Data.ZMod.Basic
import Mathlib.Topology.Algebra.Group.Basic

/-!
# tempered な位相群 —— [SemiAnbd] Definition 3.1 (i) の転写

原典: S. Mochizuki, *Semi-graphs of Anabelioids* [SemiAnbd]、
物理 p.33(`0_Source/Semi-graphs of Anabelioids.pdf`、全 92 ページ、
October 2005。**400 dpi 目視確認 2026-08-15**)。

原文 (SemiAnbd p.33, Definition 3.1):
> (i) If Π may be written as an inverse limit of an inverse system of surjections
> of countable discrete topological groups, then we shall say that Π is tempered.

原文 (SemiAnbd p.33, Remark 3.1.1):
> Observe that every profinite group is tempered.

## ★この定義は Mochizuki 固有ではない —— 明示しておく

同じページに次が書かれている。

同ページに「The notion of the tempered fundamental group is introduced in
[André], §4.」とある。★この一文は**逐語ブロックに入れられない**——`é` を
pdftotext が `Andr´e` と割るため、器具の照合が 57/63 文字で落ちる(2026-08-15 実測)。
登記済みの notationRisk (papers.json の `SemiAnbd`)そのものの実例なので、
ここは引用ではなく記述として置く。400 dpi 目視では `André` と印刷されている。

すなわち **tempered という概念そのものは André のもの**で、Mochizuki が
[SemiAnbd] §3 で行うのは「これを圏論的な視点から見直す」ことである。

原文 (SemiAnbd p.33):
> In the present manuscript, however, we wish to
> study this notion from a more categorical point of view.

したがって**本ファイルが作っているのは古典的な層であり、Mochizuki 固有の対象ではない**。
Mochizuki 固有なのは temperoid(Definition 3.1 (ii))と、半グラフ的アナベリオイドの
tempered 基本群 `π₁^temp(G)`(Proposition 3.6、物理 p.38)の側である。
そちらには**届いていない**。理由は下の「届かない地点」に書く。

## 到達点

- `SurjSystem` —— 可算離散位相群の、全射からなる逆系(Definition 3.1 (i) の「逆系」)。
- `IsTempered` —— Definition 3.1 (i) の転写。
- `tateLike` / `TateLikeGroup` —— ★**非退化な証人**。
- `not_compactSpace_tateLikeGroup` —— 副有限ではない(コンパクトでない)。
- `not_discreteTopology_tateLikeGroup` —— 離散でもない。

★退化の危険は2方向ある。「自明群」「ただの副有限群」で通してしまうことと、
「ただの可算離散群」で通してしまうこと。**どちらでもない例**を作らなければ、
この定義が profinite と discrete の合併以上のものを捉えている保証がない。
`TateLikeGroup` は `ℤ`(離散、非コンパクト)と 2進整数の塔(非離散)を同時に持つので、
**コンパクトでも離散でもない**。これは Tate 曲線の tempered 基本群が
`1 → (副有限) → Π → ℤ → 1` という形を取ること —— すなわち IUT が
temperoid を必要とする理由そのもの —— の群論的な影である。
★ただし**これは Tate 曲線の基本群そのものではない**。半グラフ的アナベリオイドを
作っていないので、そこは作れていない。

## 届かない地点(mathlib 実測、2026-08-15)

- `Mathlib/CategoryTheory/Galois/` は 11 ファイルあり、`PreGaloisCategory` /
  `FiberFunctor` / `IsFundamentalGroup` が揃っている。しかし
  **`FiberFunctor` の型は `C ⥤ FintypeCat` に焼き込まれており**、
  `PreGaloisCategory` は `HasFiniteCoproducts` と「**有限**群による商」を要求する。
  temperoid は Definition 3.1 (iii) で「**countable** colimits」を保つ関手を要求し、
  `𝓑^temp(Π)` の対象は「**countable** discrete sets」である。有限→可算は
  パラメータの変更ではないので、この API は temperoid に**インスタンス化できない**。
- `IsFundamentalGroup` は **`[CompactSpace G]` を要求する**
  (`Mathlib/CategoryTheory/Galois/IsFundamentalgroup.lean`)。
  tempered 群はまさにコンパクト性を落とすための概念なので、ここで衝突する。
- `semi-graph` は mathlib に **0 件**(`grep -rniE "semi-?graph"`)。
- 「離散群の逆極限」に相当する概念も **0 件**
  (`grep -rniE "prodiscrete|pro-discrete|inverse limit of discrete"`)。
  `Mathlib/Topology/Algebra/Category/ProfiniteGrp/Limits.lean` は
  `P ≃ₜ* limit (diagram P)` を持つが、商は `FiniteGrp` に固定されている。
-/

namespace ABC3.Found.SemiAnbd

universe u

/-- 可算離散位相群の、**全射からなる逆系**。

原文 (SemiAnbd p.33, Definition 3.1 (i)):
> (i) If Π may be written as an inverse limit of an inverse system of surjections
> of countable discrete topological groups, then we shall say that Π is tempered. -/
structure SurjSystem where
  /-- 添字。 -/
  J : Type u
  [order : Preorder J]
  [dir : IsDirected J (· ≤ ·)]
  [ne : Nonempty J]
  /-- 各段の群。 -/
  G : J → Type u
  [grp : ∀ j, Group (G j)]
  [cnt : ∀ j, Countable (G j)]
  /-- 遷移準同型。 -/
  proj : ∀ ⦃i j : J⦄, j ≤ i → G i →* G j
  /-- ★`surjections` —— 原文が要求する全射性。 -/
  proj_surj : ∀ ⦃i j : J⦄ (h : j ≤ i), Function.Surjective (proj h)
  proj_self : ∀ (j : J) (h : j ≤ j) (x : G j), proj h x = x
  proj_trans : ∀ ⦃i j k : J⦄ (h₁ : j ≤ i) (h₂ : k ≤ j) (h₃ : k ≤ i) (x : G i),
      proj h₂ (proj h₁ x) = proj h₃ x

attribute [instance] SurjSystem.order SurjSystem.dir SurjSystem.ne
  SurjSystem.grp SurjSystem.cnt

/-- ★離散性は**パラメータではなく要件**(原文 "countable **discrete** topological groups")
なので、位相は構造のフィールドにせず `⊥` を焼き込む。 -/
instance SurjSystem.instTop (S : SurjSystem.{u}) (j : S.J) : TopologicalSpace (S.G j) := ⊥

instance SurjSystem.instDisc (S : SurjSystem.{u}) (j : S.J) : DiscreteTopology (S.G j) := ⟨rfl⟩

/-- 逆系の**逆極限** —— 両立する族のなす部分群。 -/
def SurjSystem.limitSubgroup (S : SurjSystem.{u}) : Subgroup (∀ j, S.G j) where
  carrier := {x | ∀ ⦃i j : S.J⦄ (h : j ≤ i), S.proj h (x i) = x j}
  one_mem' := by intro i j h; simp
  mul_mem' := by intro a b ha hb i j h; simp [ha h, hb h]
  inv_mem' := by intro a ha i j h; simp [ha h]

/-- ★**Definition 3.1 (i) の転写**: `P` が tempered であるとは、可算離散位相群の
全射からなる逆系の逆極限として書けること。 -/
def IsTempered (P : Type u) [Group P] [TopologicalSpace P] : Prop :=
  ∃ S : SurjSystem.{u}, Nonempty (P ≃ₜ* ↥S.limitSubgroup)

/-! ### 証人となる逆系 -/

/-- 遷移準同型 `ℤ × ℤ/2^m → ℤ × ℤ/2^n`(`n ≤ m`)。第1成分は恒等、第2成分は還元。 -/
def tateProj {n m : ℕ} (h : n ≤ m) :
    Multiplicative (ℤ × ZMod (2 ^ m)) →* Multiplicative (ℤ × ZMod (2 ^ n)) :=
  AddMonoidHom.toMultiplicative
    ((AddMonoidHom.id ℤ).prodMap
      (ZMod.castHom (pow_dvd_pow 2 h) (ZMod (2 ^ n))).toAddMonoidHom)

theorem toAdd_tateProj {n m : ℕ} (h : n ≤ m) (x : Multiplicative (ℤ × ZMod (2 ^ m))) :
    Multiplicative.toAdd (tateProj h x)
      = ((Multiplicative.toAdd x).1,
          ZMod.castHom (pow_dvd_pow 2 h) (ZMod (2 ^ n)) (Multiplicative.toAdd x).2) :=
  rfl

/-- 証人の各段: `ℤ × ℤ/2ⁿ`(乗法表記)。 -/
@[reducible] def tateG (n : ℕ) : Type := Multiplicative (ℤ × ZMod (2 ^ n))

instance (n : ℕ) : Countable (tateG n) := by
  haveI : NeZero ((2 : ℕ) ^ n) := ⟨by positivity⟩
  exact inferInstanceAs (Countable (ℤ × ZMod (2 ^ n)))

/-- 各段の元を作る。 -/
def mkG (n : ℕ) (x : ℤ × ZMod (2 ^ n)) : tateG n := Multiplicative.ofAdd x

/-- 各段の元をほどく。 -/
def unG (n : ℕ) (x : tateG n) : ℤ × ZMod (2 ^ n) := Multiplicative.toAdd x

theorem unG_mkG (n : ℕ) (x : ℤ × ZMod (2 ^ n)) : unG n (mkG n x) = x := rfl

theorem unG_injective (n : ℕ) : Function.Injective (unG n) := fun _ _ h => h

theorem unG_one (n : ℕ) : unG n 1 = 0 := rfl

theorem unG_tateProj {n m : ℕ} (h : n ≤ m) (x : tateG m) :
    unG n (tateProj h x)
      = ((unG m x).1, ZMod.castHom (pow_dvd_pow 2 h) (ZMod (2 ^ n)) (unG m x).2) :=
  rfl

/-- ★**証人**: `n ↦ ℤ × ℤ/2ⁿ` からなる逆系。

第1成分の `ℤ` が「コンパクトでない」を、第2成分の 2 進の塔が「離散でない」を担う。 -/
def tateLike : SurjSystem.{0} where
  J := ℕ
  G := tateG
  proj := fun _ _ h => tateProj h
  proj_surj := by
    intro i j h
    exact Function.Surjective.prodMap Function.surjective_id
      (ZMod.castHom_surjective (pow_dvd_pow 2 h))
  proj_self := by
    intro j h x
    apply unG_injective j
    rw [unG_tateProj,
      show ZMod.castHom (pow_dvd_pow 2 h) (ZMod (2 ^ j)) = RingHom.id _ from
        Subsingleton.elim _ _]
    simp
  proj_trans := by
    intro i j k h₁ h₂ h₃ x
    apply unG_injective k
    rw [unG_tateProj, unG_tateProj, unG_tateProj]
    refine Prod.ext rfl ?_
    exact DFunLike.congr_fun
      (Subsingleton.elim
        ((ZMod.castHom (pow_dvd_pow 2 h₂) (ZMod (2 ^ k))).comp
          (ZMod.castHom (pow_dvd_pow 2 h₁) (ZMod (2 ^ j))))
        (ZMod.castHom (pow_dvd_pow 2 h₃) (ZMod (2 ^ k)))) _

/-- `ℕ` と添字 `tateLike.J` の同一視。 -/
def ix (n : ℕ) : tateLike.J := n

theorem tateLike_G (n : ℕ) : tateLike.G (ix n) = tateG n := rfl

/-- 添字を `ℕ` に戻す。 -/
def nx (j : tateLike.J) : ℕ := j

/-- ★**非退化な証人となる群**。 -/
abbrev TateLikeGroup : Type := ↥tateLike.limitSubgroup

theorem isTempered_tateLikeGroup : IsTempered TateLikeGroup :=
  ⟨tateLike, ⟨ContinuousMulEquiv.refl _⟩⟩

/-! ### 非退化性

★以下、両立条件は**すべて `ℕ` 上の独立した補題として**証明し、部分群の元を作るときに
渡す。`tateLike.J` は `ℕ` と定義的に等しいが、タクティクブロックの中では
`instances` 透明度で展開されないため、証明の中で `2 ^ j` のような式が型検査に落ちる。
これは Lean の透明度の問題であって数学の問題ではない。 -/

/-- 第 `n` 成分への射影。 -/
def ev (n : ℕ) : TateLikeGroup → tateLike.G (ix n) := fun x => (x : ∀ j, tateLike.G j) (ix n)

theorem continuous_ev (n : ℕ) : Continuous (ev n) :=
  (continuous_apply (ix n)).comp continuous_subtype_val

/-- `ℤ` の埋め込みの両立条件(`ℕ` 上で述べる)。 -/
theorem emb_compat (a : ℤ) :
    ∀ ⦃i j : ℕ⦄ (h : j ≤ i), tateProj h (mkG i (a, 0)) = mkG j (a, 0) := by
  intro i j h
  apply unG_injective j
  rw [unG_tateProj, unG_mkG, unG_mkG]
  simp

/-- 埋め込みの台となる関数(添字を `ℕ` に固定する)。 -/
def embFun (a : ℤ) : ∀ n : ℕ, tateG n := fun n => mkG n (a, 0)

/-- `ℤ` を第1成分に埋め込んだ元。 -/
def emb (a : ℤ) : TateLikeGroup := ⟨embFun a, emb_compat a⟩

theorem unG_ev_emb (a : ℤ) : unG 0 (ev 0 (emb a)) = (a, 0) := rfl

/-- ★**コンパクトでない** —— したがって副有限ではない。 -/
theorem not_compactSpace_tateLikeGroup : ¬ CompactSpace TateLikeGroup := by
  intro hc
  have hfin : (Set.range (ev 0)).Finite :=
    (isCompact_range (continuous_ev 0)).finite_of_discrete
  refine (Set.infinite_of_injective_forall_mem
    (f := fun a : ℤ => ev 0 (emb a)) ?_ (fun a => Set.mem_range_self _)) hfin
  intro a b hab
  have h2 : unG 0 (ev 0 (emb a)) = unG 0 (ev 0 (emb b)) := congrArg _ hab
  rw [unG_ev_emb, unG_ev_emb] at h2
  exact (Prod.ext_iff.mp h2).1

/-- 収束列の両立条件(`ℕ` 上で述べる)。 -/
theorem seqY_compat (N : ℕ) :
    ∀ ⦃i j : ℕ⦄ (h : j ≤ i),
      tateProj h (mkG i (0, ((2 ^ N : ℕ) : ZMod (2 ^ i))))
        = mkG j (0, ((2 ^ N : ℕ) : ZMod (2 ^ j))) := by
  intro i j h
  apply unG_injective j
  rw [unG_tateProj, unG_mkG, unG_mkG]
  refine Prod.ext rfl ?_
  simpa using map_natCast (ZMod.castHom (pow_dvd_pow 2 h) (ZMod (2 ^ j))) (2 ^ N)

/-- 収束列の台となる関数(添字を `ℕ` に固定する)。 -/
def seqYFun (N : ℕ) : ∀ n : ℕ, tateG n :=
  fun n => mkG n (0, ((2 ^ N : ℕ) : ZMod (2 ^ n)))

/-- 1 に収束するが 1 でない列: 第 `n` 成分は `2^N mod 2^n`。 -/
def seqY (N : ℕ) : TateLikeGroup := ⟨seqYFun N, seqY_compat N⟩

theorem unG_ev_seqY (N n : ℕ) :
    unG n (ev n (seqY N)) = (0, ((2 ^ N : ℕ) : ZMod (2 ^ n))) := rfl

theorem seqY_ne_one (N : ℕ) : seqY N ≠ 1 := by
  intro hN
  have h1 : unG (N + 1) (ev (N + 1) (seqY N)) = unG (N + 1) (ev (N + 1) 1) :=
    congrArg _ (congrArg _ hN)
  rw [unG_ev_seqY] at h1
  have h0 : ((2 ^ N : ℕ) : ZMod (2 ^ (N + 1))) = 0 := (Prod.ext_iff.mp h1).2
  rw [(ZMod.charP (2 ^ (N + 1))).cast_eq_zero_iff] at h0
  have h2 : N + 1 ≤ N := (Nat.pow_dvd_pow_iff_le_right (by norm_num)).mp h0
  omega

/-- 各成分では、`N` が大きければ `seqY N` は 1 に等しい。 -/
theorem seqY_coord_eq_one {n N : ℕ} (hN : n ≤ N) : ev n (seqY N) = ev n 1 := by
  apply unG_injective n
  rw [unG_ev_seqY]
  have h0 : ((2 ^ N : ℕ) : ZMod (2 ^ n)) = 0 := by
    rw [(ZMod.charP (2 ^ n)).cast_eq_zero_iff]
    exact pow_dvd_pow 2 hN
  rw [h0]
  rfl

theorem tendsto_seqY : Filter.Tendsto seqY Filter.atTop (nhds 1) := by
  rw [Topology.IsInducing.subtypeVal.tendsto_nhds_iff, tendsto_pi_nhds]
  intro j
  rw [nhds_discrete, Filter.tendsto_pure]
  filter_upwards [Filter.eventually_ge_atTop (nx j)] with N hN
  exact seqY_coord_eq_one hN

/-- ★**離散でない**。 -/
theorem not_discreteTopology_tateLikeGroup : ¬ DiscreteTopology TateLikeGroup := by
  intro hd
  haveI : DiscreteTopology TateLikeGroup := hd
  have h := tendsto_seqY
  rw [nhds_discrete, Filter.tendsto_pure] at h
  obtain ⟨N, hN⟩ := h.exists
  exact seqY_ne_one N hN

end ABC3.Found.SemiAnbd

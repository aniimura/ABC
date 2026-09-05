import ABC3.Skeleton.PGC.Section4
import ABC3.Found.PGC.UnramifiedSubextension
import ABC3.Found.PGC.QpResidueField

/-!
# [pGC] Theorem 4.2 の現在の形は**偽**——`Φ` が自由な関数だから

`Skeleton/PGC/Section4.lean::theorem_4_2` は、原文

> the natural morphism Isom_{Q_p}(K, K′) → Out_Filt(Γ_K, Γ_K′) ... is a bijection

の「自然な射」を構成する代わりに、**パラメータ `Φ` として受け取って**
その全単射性を主張していた:

```
theorem theorem_4_2 (RF) (K K')
    (Φ : (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) → FilteredGroup.OuterIso (RF.filt K) (RF.filt K')) :
    Function.Bijective Φ
```

`Φ` に**自然性の制約が一切無い**ので、これは任意の関数について全単射を
主張していることになる。**定数関数**を1つ与えれば落ちる:

* `RF` は退化させてよい(`Gv ≡ ⊤`、`degenerateRF`)。
* `K` として `ℚ_p` の**不分岐2次拡大**を取る
  (`exists_isCyclic_gal (selfField p) 2` ——本日構築した不分岐拡大理論で
  存在と `Nat.card Gal = 2` が出る)。すると
  `Isom_{Q_p}(K,K) = Gal(K/ℚ_p)` は 2 元。
* `Φ := fun _ => 恒等外部同型` は単射でない。

したがって `theorem_4_2` は **`sorry` が埋まらないのではなく、偽**。

## 何を直すべきか

原文が「**the natural** morphism」と言っている以上、`Φ` は構成されねば
ならない。部品は既にある:

* `Found/PGC/GaloisTransfer.lean::galMulEquivOf` ——`α : K ≃ₐ[ℚ_p] K′` を
  代数閉包へ延長して共役で `Γ_K ≃* Γ_K′` を得る構成
* `galMulEquivOf_indep` ——延長の取り方によらず内部自己同型で繋がる
  (だから**外部**同型としては一意)

残る穴は `map_Gv`(フィルトレーションを保つこと)で、これは
`Interface.PGC.RamificationFiltration` に「体の同型から誘導される共役が
`Gv` を保つ」という**自然性の公理が無い**ことに帰着する
(`memory/pgc-ramification-naturality-gap.md`、`Section4.lean` の
`.needs` にも記録済み)。本物の高次分岐群なら成り立つ性質だが、
現在の `Interface` では課されていない。

★この検査は `Check/PGC/InertiaDegeneracy.lean`(`I_K` を自由データに
置くと退化する)と同じ型の発見である——**自由なデータ/自由な射は、
主張を偽にするか自明にするかのどちらかになる**。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC
open scoped NormedField Valued

/-- 退化した高次分岐フィルトレーション(すべて `⊤`)。`Interface` の3条件
(閉・正規・反単調)はすべて自明に満たされる。 -/
noncomputable def degenerateRF (p : ℕ) [Fact p.Prime] : RamificationFiltration p where
  Gv _ _ := ⊤
  isClosed _ _ := by simp
  isNormal _ _ := inferInstance
  antitone _ := fun _ _ _ => le_rfl

/-- 恒等同型はフィルトレーションを保つ。 -/
def idFilteredIso (A : FilteredGroup) : FilteredGroup.Iso A A where
  equiv := ContinuousMulEquiv.refl A.G
  map_Gv v := by
    ext g
    constructor
    · rintro ⟨h, hh, rfl⟩; exact hh
    · intro hg; exact ⟨g, hg, rfl⟩

/-- 恒等外部同型——`OuterIso A A` は空でない。 -/
def idOuterIso (A : FilteredGroup) : FilteredGroup.OuterIso A A :=
  Quotient.mk _ (idFilteredIso A)

/-- **`ℚ_p` 上の自己同型が 2 つある p進局所体が存在する**——`ℚ_p` の
不分岐2次拡大。`Found/PGC/UnramifiedExtension.lean::exists_isCyclic_gal`
(各次数に不分岐拡大がちょうど一つ、Galois 群は巡回で位数はその次数)と
`Found/PGC/AdjoinFieldConstruction.lean::adjoinField`(中間体を
`PAdicLocalField` として見る)から。 -/
theorem exists_nontrivial_gal (p : ℕ) [Fact p.Prime] :
    ∃ K : PAdicLocalField p, Nontrivial (K.carrier ≃ₐ[ℚ_[p]] K.carrier) := by
  obtain ⟨x, hrank, hu, hcyc, hcard⟩ := exists_isCyclic_gal (selfField p) 2 two_ne_zero
  refine ⟨adjoinField (selfField p) x, ?_⟩
  have h2 : Nat.card ((adjoinField (selfField p) x).carrier
      ≃ₐ[ℚ_[p]] (adjoinField (selfField p) x).carrier) = 2 := hcard
  haveI : Finite ((adjoinField (selfField p) x).carrier
      ≃ₐ[ℚ_[p]] (adjoinField (selfField p) x).carrier) :=
    Nat.finite_of_card_ne_zero (by rw [h2]; norm_num)
  rw [← Finite.one_lt_card_iff_nontrivial, h2]
  norm_num

/-- **★★★★★[pGC] Theorem 4.2 の現在の形は偽**。

`Φ` に自然性の制約が無いので、2 元の定義域上の**定数関数**が反例になる。
`sorry` が埋まらないのではなく、埋めようがない。 -/
theorem theorem_4_2_statement_false (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (RF : RamificationFiltration p) (K K' : PAdicLocalField p)
        (Φ : (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) → FilteredGroup.OuterIso (RF.filt K) (RF.filt K')),
      Function.Bijective Φ) := by
  intro h
  obtain ⟨K, hK⟩ := exists_nontrivial_gal p
  obtain ⟨a, b, hab⟩ := hK.exists_pair_ne
  exact hab ((h (degenerateRF p) K K (fun _ => idOuterIso ((degenerateRF p).filt K))).1 rfl)

end ABC3.Check.PGC

import ABC3.Found.PGC.GaloisTransferContinuous
import ABC3.Found.PGC.FilteredGroup
import ABC3.Interface.PGC.LocalFieldData

/-!
# 高次分岐群の**自然性**と、Theorem 4.2 の「自然な射」

`Check/PGC/Theorem42Degenerate.lean` で示したとおり、
`Skeleton/PGC/Section4.lean::theorem_4_2` が「自然な射」を
**自由なパラメータ `Φ`** で受け取っていた形は**偽**だった
(定数関数が反例)。原文が「the **natural** morphism」と言う以上、
射は構成されねばならない。

## 部品の在庫(すべて構築済み)

| 部品 | 供給元 |
|---|---|
| 延長 `ᾱ : K̄ ≃+* K̄′` | `GaloisTransfer.lean::extendToClosure`(mathlib `IsAlgClosure.equivOfEquiv`) |
| 共役 `Γ_K ≃* Γ_{K′}` | `GaloisTransfer.lean::galMulEquivOf` |
| 延長の選択に依らない(外部同型としては一意) | `GaloisTransfer.lean::galMulEquivOf_indep` |
| 連続性 | `GaloisTransferContinuous.lean::galContinuousMulEquiv` |
| **`map_Gv`(フィルトレーション保存)** | ← **無い** |

## 残る穴は `Interface` の側にある

`Interface.PGC.RamificationFiltration` は `K` ごとに完全に独立な `Gv` を
許す抽象データであり、「体の同型から誘導される共役が `Gv` を保つ」という
**自然性の公理を持たない**(`memory/pgc-ramification-naturality-gap.md`)。
本物の高次分岐群(上付き番号付け)なら成り立つ性質だが、Herbrand の定理が
mathlib に無いため `RamificationFiltration` はまだ抽象データのままである。

そこで本ファイルは自然性を **`IsNaturalFiltration` として
明示的に切り出し**、それを仮定したうえで自然な射
`Isom_{Q_p}(K,K′) → Out_Filt(Γ_K, Γ_{K′})` を**構成する**。
`Skeleton/PGC/Section4.lean::theorem_4_2` はこの構成された射について
述べるように直した——これで「任意の関数が全単射」という偽の主張ではなく、
原文どおりの主張になる。

★非空虚性: 退化したフィルトレーション(`Gv ≡ ⊤`)は `IsNaturalFiltration` を満たす
(`exists_isNaturalFiltration`)ので、この仮説は空虚ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `RamificationFiltration` から `FilteredGroup` を作る。

★`Skeleton/PGC/Section3.lean::filteredGroupOf` および
`Skeleton/PGC/Section4.lean::RamificationFiltration.filt` と**同じ構造体
リテラル**なので定義的に等しい。ここに独立に置くのは import の向き
(`Found` は `Skeleton/PGC/Section3` 以降を import しない)のため。 -/
noncomputable def filtOf (RF : RamificationFiltration p) (K : PAdicLocalField p) :
    FilteredGroup :=
  { G := K.absGal, Gv := RF.Gv K, isClosed := RF.isClosed K, isNormal := RF.isNormal K,
    antitone := RF.antitone K }

/-- **高次分岐群の自然性**——体の同型 `α : K ≃ₐ[ℚ_p] K′` から誘導される
`Γ_K ≃ Γ_{K′}` が、各 `v` で `Γ^v` を `Γ'^v` に写す。

本物の高次分岐群なら成り立つが、`Interface.PGC.RamificationFiltration` は
これを課していない——だから Theorem 4.2 では**明示的な仮説**にする。

★名前を `RamificationFiltration.IsNatural` にしないのは、`namespace
ABC3.Found.PGC` の内側でそう書くと `ABC3.Found.PGC.RamificationFiltration.IsNatural`
になってしまうため(`Skeleton/PGC/Section4.lean` 冒頭の注意と同じ罠)。 -/
def IsNaturalFiltration (RF : RamificationFiltration p) : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (v : ℝ),
    Subgroup.map (galContinuousMulEquiv α).toMulEquiv (RF.Gv K v) = RF.Gv K' v

/-- 退化したフィルトレーション `Gv ≡ ⊤` は自然性を満たす——`IsNatural` は
空虚な仮説ではない(G2 の非空虚性対照)。 -/
noncomputable def topFiltration (p : ℕ) [Fact p.Prime] : RamificationFiltration p where
  Gv _ _ := ⊤
  isClosed _ _ := by simp
  isNormal _ _ := inferInstance
  antitone _ := fun _ _ _ => le_rfl

theorem exists_isNaturalFiltration (p : ℕ) [Fact p.Prime] :
    ∃ RF : RamificationFiltration p, IsNaturalFiltration RF := by
  refine ⟨topFiltration p, ?_⟩
  intro K K' α v
  exact Subgroup.map_top_of_surjective _ (galContinuousMulEquiv α).surjective

/-- **★★★★★自然な射(filtered group の同型版)**——体の同型 `α` から、
高次分岐群のフィルトレーションを保つ `Γ_K ≅ Γ_{K′}` を作る。 -/
noncomputable def naturalFilteredIso (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
    {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.Iso (filtOf RF K) (filtOf RF K') where
  equiv := galContinuousMulEquiv α
  map_Gv v := hnat α v

/-- **★★★★★★[pGC] Theorem 4.2 の「自然な射」**
`Isom_{Q_p}(K, K′) → Out_Filt(Γ_K, Γ_{K′})`。

`galMulEquivOf_indep`(延長の取り方は内部自己同型のずれしか生まない)が
あるので、**外部**同型としては延長の選択に依らない。 -/
noncomputable def naturalOuterIso (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
    {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.OuterIso (filtOf RF K) (filtOf RF K') :=
  Quotient.mk _ (naturalFilteredIso RF hnat α)

end ABC3.Found.PGC

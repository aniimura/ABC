import ABC3.Found.PGC.LubinTateTowerCompatible

/-!
# 無限に続く compatible な原始的生成元の列 `(x_n)_{n≥1}`(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の最終組み立てに要る「無限
compatible列」——`x_n:=π·x_{n+1}`(全`n`で成り立つ)を満たす、`ψ_n`
の根(原始的なπ^n-捩れ点)の無限列——を、`Found/PGC/LubinTateTower
Compatible.lean::exists_pi_mul_eq_of_mem_iteratedLubinTatePsiTorsion
Points`(「1段分の全射性」、既出)を`Nat.rec`で積み重ねて実際に
構成する。

## `PsiGen`(状態の束ね方)

`level m+1`の`ψ_{m+1}`の根`x`と、それに付随する諸事実(`Λ_{m+1}`
所属・`x∈K.carrier⟮x⟯`・有限次元性)を1つの構造`PsiGen K...m`に
束ねる(`m`は0始まり、`level:=m+1`)。

## 帰納法の1ステップと、そこで踏んだ罠

`psiGenStepEx`→`PsiGenStepResult`→`psiGenStepResult`→`psiGenStep`
という順に組む。★★★ここで**新しい罠**を踏んだ:「独立に書いた型
注釈と、既存の項の型を照合させる」ことが、`Classical.choice`由来の
深いnested `Exists.choose`を経由する項に対しては**極端に高くつく**
(`maxHeartbeats`を100万まで上げても`whnf`のtimeoutが取れなかった)。
個々の射影(`.pt`だけの等式等)は`rfl`で一瞬(0.3秒)で通るのに、
複数の射影を組み合わせた**独立な型注釈**を書いて`exact`で照合させる
と刺さる、という非対称な挙動を実測した。

**回避策**: 型注釈を省略し(`theorem`ではなく`noncomputable def`を
使い、右辺の項からLeanに型を**推論させる**)、Leanに「照合」させる
のではなく「読み取らせる」ことで、この罠を完全に回避した——
`psiGenStep_compat`・`psiGenSeq_compat`はどちらもこの形。得られる
型は(見た目は少し遠回りだが)`psiGenStepResult`の`pt`/`hn`/`hmem`を
経由する形になるが、これらは個々に`rfl`で`psiGenStep`/`psiGenSeq`の
対応する射影と一致することを確認済み——使う側は必要な瞬間に
`rfl`で橋渡しすればよい。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- level `m+1` の `ψ_{m+1}` の根(原始的なπ^{m+1}-捩れ点)と、それに
付随する諸事実を束ねたもの。 -/
structure PsiGen {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) where
  pt : K.closure
  hψ : pt ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
  hn : pt ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
  hmem : pt ∈ IntermediateField.adjoin K.carrier ({pt} : Set K.closure)
  hfd : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({pt} : Set K.closure))

/-- 出発点(`m=0`、level `1`): `ψ_1` の根の非空性(既出)から任意に
1つ選ぶ。 -/
noncomputable def psiGenBase {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff)) :
    PsiGen K hq hπmax hπne0 f hf0 hf1 hf 0 :=
  have hne := iteratedLubinTatePsiTorsionPoints_nonempty K hq hπmax hπne0 f hf0 hf1 hf 1 (le_refl 1)
  have hψ : hne.choose ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (0 + 1) (by omega) :=
    hne.choose_spec
  have hn : hne.choose ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (0 + 1) := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf 1 (le_refl 1)]
    exact Finset.mem_union_right _ hψ
  have hmem : hne.choose ∈ IntermediateField.adjoin K.carrier ({hne.choose} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier hne.choose
  { pt := hne.choose
    hψ := hψ
    hn := hn
    hmem := hmem
    hfd := finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf
      (0 + 1) hne.choose hn }

/-- 帰納法の1ステップの出力——「次の生成元」の諸事実と、`π·(次の
生成元)=(前の生成元)`という compat 条件を**1つの構造にまとめて**
束ねる(`hcompat`が`pt`/`hn`/`hmem`という**兄弟フィールド**を参照
する形——これが罠を避ける鍵。`gen.pt`のような**入れ子越しの**参照
にすると`whnf`のtimeoutを踏むことを実測した)。 -/
structure PsiGenStepResult {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (prevPt : K.closure) where
  pt : K.closure
  hψ : pt ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) (by omega)
  hn : pt ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1)
  hmem : pt ∈ IntermediateField.adjoin K.carrier ({pt} : Set K.closure)
  hfd : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({pt} : Set K.closure))
  hcompat : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) pt hn hmem π) :
      IntermediateField.adjoin K.carrier ({pt} : Set K.closure)) : K.closure) = prevPt

/-- 帰納法の1ステップの本体: `exists_pi_mul_eq_of_mem_iteratedLubinTate
PsiTorsionPoints`(既出)を、新しい種(`ψ_{m+2}`の非空性から)と
`prev.pt`(目標`y`)を使って適用するだけ。 -/
noncomputable def psiGenStepResult {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (prev : PsiGen K hq hπmax hπne0 f hf0 hf1 hf m) :
    PsiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m prev.pt := by
  have hne := iteratedLubinTatePsiTorsionPoints_nonempty K hq hπmax hπne0 f hf0 hf1 hf (m + 2) (by omega)
  have hseedψ : hne.choose ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf ((m + 1) + 1)
      (by omega) := hne.choose_spec
  have hseedn : hne.choose ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf ((m + 1) + 1) := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf (m + 2) (by omega)]
    exact Finset.mem_union_right _ hseedψ
  have hseedmem : hne.choose ∈ IntermediateField.adjoin K.carrier ({hne.choose} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier hne.choose
  haveI hseedfd : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({hne.choose} : Set K.closure)) :=
    finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf
      ((m + 1) + 1) hne.choose hseedn
  have hex := exists_pi_mul_eq_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf
    (m + 1) (by omega) hne.choose hseedψ hseedn hseedmem prev.pt prev.hψ
  exact
    { pt := hex.choose
      hψ := hex.choose_spec.choose
      hn := hex.choose_spec.choose_spec.choose
      hmem := hex.choose_spec.choose_spec.choose_spec.choose
      hfd := hex.choose_spec.choose_spec.choose_spec.choose_spec.choose
      hcompat := hex.choose_spec.choose_spec.choose_spec.choose_spec.choose_spec }

/-- 「次の生成元」を`PsiGen`の形へ落とす——`PsiGenStepResult`から
`hcompat`を捨てて`pt`/`hψ`/`hn`/`hmem`/`hfd`だけを取り出すだけ。 -/
noncomputable def psiGenStep {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (prev : PsiGen K hq hπmax hπne0 f hf0 hf1 hf m) :
    PsiGen K hq hπmax hπne0 f hf0 hf1 hf (m + 1) :=
  let r := psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m prev
  { pt := r.pt, hψ := r.hψ, hn := r.hn, hmem := r.hmem, hfd := r.hfd }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**帰納法の1ステップのcompat性**: `π·(次の生成元)=(前の生成元)`。
型注釈を省略し(`noncomputable def`)、`PsiGenStepResult.hcompat`から
直接推論させることで、独立な型注釈同士を照合させる際の`whnf`
timeoutの罠を回避した。 -/
noncomputable def psiGenStep_compat {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (prev : PsiGen K hq hπmax hπne0 f hf0 hf1 hf m) :=
  (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m prev).hcompat

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**無限に続く compatible な原始的生成元の列そのもの**: 各`m`について
`level m+1`の`ψ_{m+1}`の根`(psiGenSeq m).pt`を与える。節目(5)
(射影極限)の核心的な構成物。 -/
noncomputable def psiGenSeq {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff)) :
    (m : ℕ) → PsiGen K hq hπmax hπne0 f hf0 hf1 hf m
  | 0 => psiGenBase K hq hπmax hπne0 f hf0 hf1 hf
  | (m + 1) => psiGenStep K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**列全体のcompat性**: `π·(psiGenSeq (m+1)への値)=(psiGenSeq m).pt`。
`psiGenStep_compat`と同じ理由で型注釈を省略している——得られる型は
`psiGenStepResult K...m (psiGenSeq m)`の`pt`/`hn`/`hmem`を経由する
形になるが、これらは(`psiGenStep`/`psiGenSeq`の定義から)個々に
`rfl`で`(psiGenSeq (m+1))`の対応するフィールドと一致することが
確認できる(例: `(psiGenSeq K...(m+1)).pt = (psiGenStepResult K...m
(psiGenSeq K...m)).pt`は`rfl`で通る)——使う側は必要な瞬間にこの
`rfl`で橋渡しすればよい。 -/
noncomputable def psiGenSeq_compat {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) :=
  psiGenStep_compat K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)

end ABC3.Found.PGC

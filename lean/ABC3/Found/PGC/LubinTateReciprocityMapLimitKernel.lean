import ABC3.Found.PGC.LubinTateReciprocityMapLimitSurjective

/-!
# `reciprocityMapLimit` の核の特徴づけ(`sorry` 無し)——「単射性」の
正しい言い換え

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の部品(iii)の後半。軌道修正の
記録(`memory/pgc-lubin-tate-existence-progress.md`)で確認した
とおり、`reciprocityMapLimit:Gal(K.closure/K.carrier)→*Compatible
Units`は定義域`Gal(K.closure/K.carrier)`(完全な絶対Galois群、
非可換)から値域(アーベル群)への準同型なので**単射にはなり得ない**
——正しい形式化対象は**核の特徴づけ**である。

## 結論

`σ∈ker(reciprocityMapLimit) ↔ ∀m,σ((psiGenSeq m).pt)=(psiGenSeq
m).pt`——`σ`が無限compatible列`(psiGenSeq m).pt`(古典的なLubin-Tate
理論の`Λ_∞`の生成元の列)を**すべて固定する**ことと同値。これは
古典論の`ker=Gal(K̄/K_π)`(`K_π:=K(Λ_∞)`)を、`K_π`という中間体
オブジェクトを構成せずに、その**生成元の言葉で直接**述べたもの。

## 証明の構造

(→) `reciprocityMapLimit σ=1`から`reciprocityMapLimitFamily σ(m+1)
=1`(各m)を取り出し、`principalUnitsQuotientEquiv`の単射性
(`MulEquiv.map_eq_one_iff`)で`reciprocityMap(x_m)(σ)=1`(`(𝒪_K)^×
⧸principalUnits`の単位元)、`QuotientGroup.mk_one`で`=QuotientGroup.
mk 1`に書き換えてから`reciprocityMap_spec`+`unitActionQuotientLift_
mk`+`lubinTateActionAtTorsionPoint_one`(既出、単位元の作用は恒等)
で`σ(x_m)=x_m`を得る。

(←) 逆に`∀m,σ(x_m)=x_m`から、`lubinTateActionAtTorsionPoint_one`で
`x_m=1·x_m`と書き換え、`reciprocityMap_eq_mk_of_apply_eq`(既出、
前回)で`reciprocityMap(x_m)(σ)=QuotientGroup.mk 1=1`、`principal
UnitsQuotientEquiv`の`map_one`で`reciprocityMapLimitFamily σ(m+1)
=1`。`n=0`成分は自明群なので自動的に一致。

新しい数学的な難所は無かった——既に確立済みの部品(全射性の証明で
使ったのと同じ道具立て)を組み合わせるだけで届いた。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMapLimit`の核の特徴づけ**——単射性に代わる正しい主張。
`σ`が核に入ることと、`σ`が無限compatible生成元列`(psiGenSeq m).pt`
をすべて固定することが同値。 -/
theorem mem_ker_reciprocityMapLimit_iff
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    σ ∈ MonoidHom.ker (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf) ↔
      ∀ m : ℕ, σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt := by
  haveI : ∀ m, FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) :=
    fun m => (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  have h0sub : Subsingleton (𝒪[K.carrier] ⧸ Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])))ˣ := by
    have h00 : Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])) = ⊤ := by
      rw [pow_zero]; exact Ideal.span_singleton_one
    rw [h00]
    constructor
    intro a b
    apply Units.ext
    exact Subsingleton.elim _ _
  rw [MonoidHom.mem_ker]
  constructor
  · intro hσ m
    have h1 : reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1) = 1 :=
      congrFun (congrArg Subtype.val hσ) (m + 1)
    change principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) = 1 at h1
    rw [MulEquiv.map_eq_one_iff] at h1
    have hspec := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
    rw [h1, ← QuotientGroup.mk_one, unitActionQuotientLift_mk] at hspec
    rw [← hspec, show ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) = 1 from rfl,
      lubinTateActionAtTorsionPoint_one]
  · intro hσ
    apply Subtype.ext
    funext n
    match n with
    | 0 => exact h0sub.elim _ _
    | (m + 1) =>
      show principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
          (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) = 1
      have hσm := hσ m
      have h1 : σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
          (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])) :
            IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) : K.closure) := by
        rw [show ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) = 1 from rfl, lubinTateActionAtTorsionPoint_one]
        exact hσm
      rw [reciprocityMap_eq_mk_of_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) _
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ _ _ 1 σ h1,
        QuotientGroup.mk_one]
      exact map_one _

end ABC3.Found.PGC

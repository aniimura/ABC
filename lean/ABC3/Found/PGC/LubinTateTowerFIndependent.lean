import ABC3.Found.PGC.DworkThetaEval

/-!
# `K̂^m_f = K̂^m_{f′}` —— 添字 `f` が落ちる最初の一歩(Yoshida 2008 Corollary 4.9 前半)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Corollary 4.9(物理 p.9)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#cor-4-9`。

原文 (Yoshida08 p.9):
> Corollary 4.9. The K^m_f and ρ_f,m, hence also K^LT_f and ρ_f, of Proposition 4.7(ii) do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

原典の証明(p.9、`0_Source` の `.txt` の 516–522 行を直読した全文):

> Proof. For f, f′ with linear coefficients π, π′, take θ ∈ Θ^{K̂,×}_{π,π′} and
> [θ] : µ_{f,m} ≅→ µ_{f′,m} by Proposition 4.8. Lemma 4.6 shows K̂^m_f = K̂^m_{f′}.
> If σ(α) = [xπ_j](α) for σ ∈ W(K̂^m_f/K), then σ([θ](α)) = [θ]^{(j)}[xπ_j](α) = [xπ′_j][θ](α)
> by Lemma 4.5, hence ρ_{f,m} = ρ_{f′,m}. □

★★本ノードは Corollary 4.9 の**前半のみ** —— 「Lemma 4.6 shows K̂^m_f = K̂^m_{f′}」
(体の一致)である。後半 `ρ_{f,m} = ρ_{f′,m}`(写像の一致)は **Lemma 4.5 が木に
まだ無い**ため、別ノードとして立てる(本ファイルには一切含まれない)。

原典が引く Lemma 4.6 の証明はこう書かれている(p.8):

> As [θ], [θ^{-1}] ∈ O_L[[X]], we have µ_{f′,m} = [θ](µ_{f,m}) ⊂ L^m_f and
> µ_{f,m} ⊂ L^m_{f′}, thus L^m_f = L^m_{f′}.

すなわち「`θ` が `𝒪_L` 係数の**冪級数**であるから、その値は `L^m_f` に入る」。
★ここが原典の畳み方であり、本ファイルの本体でもある: `θ(λ)` は多項式の値ではなく
**冪級数の値=極限**なので、`L^m_f = K̂^{ur}(Λ_{f,m})` が `ℂ_K` の中で**閉**である
ことが要る。これを実機で確かめた(下の「#8 `SubfieldClosed` は 0 行ではない」)。

## 何を足したか

### 抽象核(分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

| 宣言 | 内容 |
|---|---|
| `aeval_mem_of_isClosed` | ★★**閉部分環は冪級数の値を含む**: `T` が閉部分環で係数環の像と評価点を含むなら `θ(z) ∈ T`。位相環の一般論だけ(トランケーション `trunc N θ` の極限) |
| `constantCoeff_of_subst_eq_X` | `θ ∘ θ′ = X` と `θ(0) = 0` なら `θ′(0) = 0`(形式冪級数だけ) |
| `isClosed_intermediateField_of_finiteDimensional` | ★★**完備体上有限次元の中間体は閉**。`Submodule.closed_of_finiteDimensional` の中間体版 |

★`aeval_mem_of_isClosed` は「集合の間に全単射があるだけでは環は一致しない」という
問題への答えである。**要るのは全単射ではなく「その写像が θ の評価で与えられている」
という情報**で、閉部分環はその情報を吸収できる。

### 具体層

| 宣言 | 内容 |
|---|---|
| `algebra_unramifiedCompletion_closureCompletion` | `ℂ_K` を `K̂^{ur}`-代数として見る(`ι_ur` から) |
| `isIntegral_closureCompletionCoe` | `K^al` の元の `ℂ_K` での像は `K̂^{ur}` 上整 |
| `torsionCompletionSet` | `Λ_{f,n}` の `𝒪_{ℂ_K}` の中での像 |
| `torsionCompletionSet_subset_of_isClosed` | ★★**輸送核**: 閉部分環が `Λ_{g,n}` の像を含めば `Λ_{f,n}` の像も含む |
| `lubinTateCompletionField` | ★**`K̂^m_f`** = `K̂^{ur}` に `Λ_{f,n}` を添加した `ℂ_K` の中間体 |
| `finiteDimensional_lubinTateCompletionField` / `isClosed_lubinTateCompletionField` | `K̂^m_f` は `K̂^{ur}` 上有限次元、したがって `ℂ_K` の中で閉 |
| `lubinTateCompletionField_le_of_dworkTheta` | 片側の包含 `K̂^m_f ≤ K̂^m_g` |
| `lubinTateCompletionField_eq_of_dworkTheta` | ★★**両側** —— `θ` を仮定に置いた形 |
| `lubinTateCompletionField_eq` | ★★★**完成形** —— `θ` を Λ6(Dwork)から取ってくる形 |
| `completionAdjoin` 以下 | 整数側の影(`𝒪_{K̂^{ur}}[Λ_{f,n}]` の閉包)。同じ輸送核から出る |

## 「`θ(λ)` が adjoin に入る」をどう解決したか —— (a) と (b) の合わせ技

段取り係が挙げた 4 択のうち **(a) 在庫にあった + (b) 自前で書いた**である。

* `Submodule.closed_of_finiteDimensional`(mathlib、`Analysis/Normed/Module/FiniteDimension`)は
  `NormedSpace` を要求せず **`Module` + `ContinuousSMul` + `T2Space` + `IsTopologicalAddGroup`**
  だけで済む。★これが決定打だった(`NormedSpace K̂^{ur} ℂ_K` を組む必要が無い)。
* `NontriviallyNormedField (unramifiedCompletion K)` は**木にすでにある**
  (`nontriviallyNormedField_unramifiedCompletion`)。
* `IntermediateField.finiteDimensional_adjoin`(mathlib)も有限集合の整元添加でそのまま使える。

★★**したがって #8 `SubfieldClosed` の「0 行(削除提案)」は誤りである。**
「完備体の有限次拡大は完備だから自動で出る」という推測は**方向としては正しい**が、
実機では次の 5 つを自分で置く必要があった(本ファイルの §0 末尾と §1、実測 **約 45 行**):

1. `Algebra (unramifiedCompletion K) (closureCompletion K)`(`ι_ur` からの scoped instance)
2. `IsScalarTower K.carrier (unramifiedCompletion K) (closureCompletion K)`
3. `ContinuousSMul (unramifiedCompletion K) (closureCompletion K)`
4. `closureCompletionAlgHom`(`K^al →ₐ[K] ℂ_K`)と `isIntegral_closureCompletionCoe`
5. `isClosed_intermediateField_of_finiteDimensional`(抽象核、8 行)

★「0 行」ではないが「250 行」でもない。**45 行**である。

## 両方の包含が入ったか —— 入った

`lubinTateCompletionField_eq_of_dworkTheta` は `le_antisymm` で
`K̂^m_f ≤ K̂^m_g` と `K̂^m_g ≤ K̂^m_f` の**両方**を出している。逆向きは
抽象核 `subst_intertwine_of_comp_inverse`(#7 が用意した)で `θ′` の絡みを作り、
同じ補題に `(f, g, θ, θ′) := (g, f, θ′, θ)` を代入するだけである。
★`θ′` の定数項が `0` であることは Λ6 の結論に**含まれていない**ので、
`constantCoeff_of_subst_eq_X` を新たに用意して `θ ∘ θ′ = X` から導いた。

## 退化の自己検査

* ★**`f` と `g` の線形係数が `π` と `uπ`(`u` 単元)で結ばれていることを落とすと偽。**
  無関係な 2 つの Lubin-Tate 塔は一致しない。完成形
  `lubinTateCompletionField_eq` では `hg1 : coeff 1 g = u * π` と `hu : IsUnit u`
  がそれである(`θ` の存在 = Prop 4.8 = Λ6 がそもそもこの形でしか出ない)。
* ★**`θ` の存在を落とすと空虚。** `lubinTateCompletionField_eq_of_dworkTheta` は
  `θ`・`θ′`・絡み `hint` を仮定に取り、完成形はそれを
  `exists_arithFrobenius_isCoherent_dworkThetaStep2` から供給している。
* ★**`n` を固定しないと `Λ` が定まらない。** すべての宣言が `n` を明示に取る。
* ★**片側だけでは「一致」ではない。** 上記のとおり両方出している。
* ★**`K̂^m_f` が閉であることを落とすと `θ(λ) ∈ K̂^m_f` が言えない。**
  `isClosed_lubinTateCompletionField` を落とすと
  `lubinTateCompletionField_le_of_dworkTheta` の証明は
  `torsionCompletionSet_subset_of_isClosed` の `hT` を埋められない。
* ★**`aeval_mem_of_isClosed` から `IsClosed` を落とすと偽。**
  `𝒪_{K̂^{ur}}[λ]` は一般に `θ(λ)` を含まない(冪級数の値は多項式の値ではない)。

## 逸脱(記録)

1. **`K̂^m_f` の定義**。原典は `L^m_f := L(µ_{f,m})`(`L = K̂`)と書く。本ファイルは
   `ℂ_K`(`K^al` の完備化、M3)の中で
   `IntermediateField.adjoin (unramifiedCompletion K) (ι_al '' Λ_{f,n})` と定義する。
   `ι_al` は単射(`closureCompletionCoe_injective`)なので情報は落ちていない。
   ★原典の `µ_{f,m}` は `[π^m]_f` の零点集合、本ファイルの `Λ_{f,n}` は在庫の
   `iteratedLubinTateTorsionPoints`(Weierstrass 分解の `D_n` の根集合)であり、
   両者の一致は #7(`DworkThetaEval.lean`)で確認済み。
2. **添字**。原典の `m`(レベル)を本ファイルは `n` と書く(在庫の
   `iteratedLubinTateTorsionPoints` の引数名に合わせた)。
3. **整数側の影を併記した。** 原典には無いが、同じ輸送核から
   `𝒪_{K̂^{ur}}` と `Λ_{f,n}` が生成する `𝒪_{ℂ_K}` の**閉部分環**の一致
   (`completionAdjoin_torsionCompletionSet_eq`)も出しておいた。
   ★こちらは位相閉包を取った形で述べている(`𝒪_{K̂^{ur}}[Λ]` 自身が閉であること
   ——完備離散付値環上の束の閉性——は本ファイルでは扱っていない)。
   体の側 `lubinTateCompletionField_eq` には位相閉包は**出てこない**(有限次元
   だから自動的に閉)。原典に対応するのは体の側である。
4. **`scoped instance` を 2 つ足した**(`algebra_unramifiedCompletion_closureCompletion`・
   `isScalarTower_carrier_unramifiedCompletion_closureCompletion`)。どちらも
   `scoped` なので `namespace ABC3.Found.PGC` の中でしか見えず、かつ本ファイルを
   直接 import する `Found/PGC/*.lean` は現在ゼロ(`Found.lean` の集約 import のみ)
   なので、既存の instance 探索には影響しない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

def lubinTateCompletionField_eq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

/-! ### このファイル限りのインスタンス

`𝒪_{ℂ_K}` の線形位相と連続スカラー倍を**このファイル限りの**インスタンスにする
(`DworkThetaEval.lean` と同じ理由——`PowerSeries.aeval` を statement に出すので
`haveI` では間に合わない。`tools/lean-idioms.md` #105)。 -/

attribute [local instance] isLinearTopology_closureCompletionInt continuousSMul_closureCompletionInt

/-! ## 0. 抽象核

分岐・付値・Lubin-Tate の語彙が 1 つも出てこない部分。原典の設定に依らない。 -/

section AbstractCore

open scoped PowerSeries.WithPiTopology in
/-- ★★★★★★★★★★★★★★★★**閉部分環は冪級数の値を含む**。

`T` が位相環 `S` の**閉**部分環で、係数環 `A` の像と評価点 `z` を含むなら、
`z` における任意の冪級数の値 `θ(z)` は `T` に入る。

段取り: `θ(z) = lim_N (trunc N θ)(z)`(`PowerSeries.continuous_aeval` と
`tendsto_trunc_atTop`)で、各 `(trunc N θ)(z)` は `A` の像と `z` の多項式なので
`T` の元。あとは `IsClosed.mem_of_tendsto`。

★これが本ファイルの本体である: 原典 Lemma 4.6 の
「`[θ] ∈ 𝒪_L[[X]]` だから `[θ](µ_{f,m}) ⊂ L^m_f`」がまさにこの補題。
★★`IsClosed` を落とすと**偽**(冪級数の値は多項式の値ではない)。 -/
theorem aeval_mem_of_isClosed {A S : Type*} [CommRing A] [CommRing S] [Algebra A S]
    [UniformSpace A] [UniformSpace S] [IsUniformAddGroup A] [IsTopologicalSemiring A]
    [IsUniformAddGroup S] [T2Space S] [CompleteSpace S] [IsTopologicalRing S]
    [IsLinearTopology S S] [ContinuousSMul A S]
    (T : Subring S) (hT : IsClosed (T : Set S))
    (hbase : ∀ a : A, algebraMap A S a ∈ T)
    {z : S} (hz : PowerSeries.HasEval z) (hzT : z ∈ T) (θ : PowerSeries A) :
    PowerSeries.aeval hz θ ∈ T := by
  have h3 : Filter.Tendsto (fun N => PowerSeries.aeval hz
      ((PowerSeries.trunc N θ : Polynomial A) : PowerSeries A))
      Filter.atTop (nhds (PowerSeries.aeval hz θ)) :=
    (PowerSeries.continuous_aeval hz).continuousAt.tendsto.comp
      (PowerSeries.WithPiTopology.tendsto_trunc_atTop (R := A) θ)
  refine hT.mem_of_tendsto h3 (Filter.Eventually.of_forall fun N => ?_)
  rw [PowerSeries.aeval_coe]
  induction (PowerSeries.trunc N θ) using Polynomial.induction_on with
  | C a => rw [Polynomial.aeval_C]; exact hbase a
  | add P Q hP hQ => rw [map_add]; exact T.add_mem hP hQ
  | monomial k a _ =>
      rw [map_mul, Polynomial.aeval_C, map_pow, Polynomial.aeval_X]
      exact T.mul_mem (hbase a) (T.pow_mem hzT _)

/-- 合成逆の定数項は自動的に `0`——`θ ∘ θ′ = X` と `θ(0) = 0` から。

★Λ6(`DworkThetaStep2.lean`)の結論は `θ′` の定数項を**出していない**ので、
本ファイルはここで自前で作る。`θ′ = C c + X·u` と分解して `subst` を通し、
定数項を比べるだけ。 -/
theorem constantCoeff_of_subst_eq_X {A : Type*} [CommRing A] {θ θ' : PowerSeries A}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (h : PowerSeries.subst θ θ' = PowerSeries.X) :
    PowerSeries.constantCoeff θ' = 0 := by
  have hHS : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  set c := PowerSeries.constantCoeff θ' with hc
  have h1 : PowerSeries.constantCoeff (θ' - PowerSeries.C c) = 0 := by simp [hc]
  obtain ⟨u, hu⟩ := PowerSeries.X_dvd_iff.mpr h1
  have h2 : θ' = PowerSeries.C c + PowerSeries.X * u := by rw [← hu]; ring
  rw [h2, PowerSeries.subst_add hHS, PowerSeries.subst_mul hHS, PowerSeries.subst_C,
    PowerSeries.subst_X hHS] at h
  have h3 := congrArg PowerSeries.constantCoeff h
  have h4 : PowerSeries.constantCoeff ((MvPowerSeries.C c : PowerSeries A)) = c := rfl
  simpa [hθ0, h4] using h3

/-- ★★★★★★★★★★★★**完備体上有限次元の中間体は閉**。

`Submodule.closed_of_finiteDimensional`(mathlib)の中間体版。
★mathlib のこの補題は `NormedSpace` を要求せず、`Module` + `ContinuousSMul` +
`T2Space` + `IsTopologicalAddGroup` だけで済む——これが本ファイルを可能にした。

★★これが段取り係の言う **#8 `SubfieldClosed`** の中身である。
「0 行(自動で出る)」ではなく、この 8 行 + 具体層の instance 約 37 行が要った。 -/
theorem isClosed_intermediateField_of_finiteDimensional {𝕜 L : Type*}
    [NontriviallyNormedField 𝕜] [CompleteSpace 𝕜] [Field L] [TopologicalSpace L]
    [IsTopologicalAddGroup L] [T2Space L] [Algebra 𝕜 L] [ContinuousSMul 𝕜 L]
    (E : IntermediateField 𝕜 L) [FiniteDimensional 𝕜 E] :
    IsClosed (E : Set L) := by
  haveI : FiniteDimensional 𝕜 ↥(Subalgebra.toSubmodule E.toSubalgebra) :=
    inferInstanceAs (FiniteDimensional 𝕜 E)
  have h := Submodule.closed_of_finiteDimensional (Subalgebra.toSubmodule E.toSubalgebra)
  rwa [Subalgebra.coe_toSubmodule] at h

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. `ℂ_K` を `K̂^{ur}`-代数として見る

M3(`ClosureCompletion.lean`)は `𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}` の代数構造を作ったが、
**体のレベル** `K̂^{ur} → ℂ_K` の代数構造は作っていない。原典 Lemma 4.6 の
`L^m_f = L(µ_{f,m})` は体の話なのでここで足す。★どちらも `ι_ur` から作るので
`𝒪` 側と整合している(`algebraMap_unramifiedInt_closureCompletionInt_coe`)。 -/

/-- **`ℂ_K` は `K̂^{ur}`-代数**(`ι_ur : K̂^{ur} →+* ℂ_K` から)。 -/
@[reducible] noncomputable scoped instance algebra_unramifiedCompletion_closureCompletion
    (K : PAdicLocalField p) : Algebra (unramifiedCompletion K) (closureCompletion K) :=
  (unramifiedToClosureCompletion K).toAlgebra

theorem algebraMap_unramifiedCompletion_closureCompletion_eq (K : PAdicLocalField p)
    (z : unramifiedCompletion K) :
    algebraMap (unramifiedCompletion K) (closureCompletion K) z
      = unramifiedToClosureCompletion K z := rfl

/-- スカラー塔 `K → K̂^{ur} → ℂ_K`。 -/
scoped instance isScalarTower_carrier_unramifiedCompletion_closureCompletion
    (K : PAdicLocalField p) :
    IsScalarTower K.carrier (unramifiedCompletion K) (closureCompletion K) :=
  IsScalarTower.of_algebraMap_eq fun a => (unramifiedToClosureCompletion_algebraMap K a).symm

/-- `K̂^{ur}` の `ℂ_K` への作用は連続(`ι_ur` が連続で、`ℂ_K` の積が連続だから)。
★`Submodule.closed_of_finiteDimensional` が要求する 4 条件のうちの 1 つ。 -/
theorem continuousSMul_unramifiedCompletion_closureCompletion (K : PAdicLocalField p) :
    ContinuousSMul (unramifiedCompletion K) (closureCompletion K) :=
  ⟨continuous_mul.comp
    (((continuous_unramifiedToClosureCompletion K).comp continuous_fst).prodMk continuous_snd)⟩

attribute [local instance] continuousSMul_unramifiedCompletion_closureCompletion

/-- `ι_al : K^al → ℂ_K` を `K`-代数準同型として見る。 -/
noncomputable def closureCompletionAlgHom (K : PAdicLocalField p) :
    K.closure →ₐ[K.carrier] closureCompletion K :=
  { closureCompletionCoe K with commutes' := fun _ => rfl }

/-- **`K^al` の元の `ℂ_K` での像は `K̂^{ur}` 上整**。

`K^al/K` が整(`Algebra.IsIntegral K.carrier K.closure`)であることを `ι_al` で送り、
スカラー塔 `K → K̂^{ur} → ℂ_K` で底を持ち上げる(`IsIntegral.tower_top`)。
★`Λ_{f,n}` の元が `K̂^{ur}` 上代数的であること——`K̂^m_f` の有限次元性の入口。 -/
theorem isIntegral_closureCompletionCoe (K : PAdicLocalField p) (lam : K.closure) :
    IsIntegral (unramifiedCompletion K) (closureCompletionCoe K lam) :=
  ((Algebra.IsIntegral.isIntegral (R := K.carrier) lam).map (closureCompletionAlgHom K)).tower_top

/-! ## 2. 捩れ点の像と輸送核

`Λ_{f,n} ⊆ K^al` の `𝒪_{ℂ_K}` の中での像を集合として置き、
「閉部分環が `Λ_{g,n}` の像を含むなら `Λ_{f,n}` の像も含む」を 1 本にまとめる。
★この 1 本が体の側にも整数の側にも効く。 -/

/-- `Λ_{f,n}` の `𝒪_{ℂ_K}` の中での像(`ι_al` の像として述べる)。 -/
def torsionCompletionSet (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) : Set ↥(closureCompletionInt K) :=
  {w | ∃ r ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
        (w : closureCompletion K) = closureCompletionCoe K r}

/-- ★★★★★★★★★★★★★★★★★★★★**輸送核**——`𝒪_{ℂ_K}` の閉部分環 `T` が
`𝒪_{K̂^{ur}}` の像と `Λ_{g,n}` の像を含むなら、`Λ_{f,n}` の像も含む。

段取り(原典 Lemma 4.6 の 1 行そのもの):
1. `λ ∈ Λ_{f,n}` に対し #7 の `exists_mem_torsionPoints_evalAt` で
   `θ(λ) = ι_al r`、`r ∈ Λ_{g,n}`。したがって `θ(λ) ∈ T`。
2. `θ′` を `θ(λ)` で評価すると `θ′(θ(λ)) = λ`
   (#7 の `coe_aeval_evalAt_comp_inverse`、`subst θ θ′ = X` から)。
3. `θ(λ) ∈ T`・`T` は閉・`T ⊇ 𝒪_{K̂^{ur}}` なので、抽象核
   `aeval_mem_of_isClosed` から `θ′(θ(λ)) ∈ T`、すなわち `ι_al λ ∈ T`。

★`f` と `g` は別々の素元 `π`・`ϖ` に属してよい(`hint` が両者を繋ぐ)。 -/
theorem torsionCompletionSet_subset_of_isClosed (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) (T : Subring ↥(closureCompletionInt K))
    (hT : IsClosed (T : Set ↥(closureCompletionInt K)))
    (hbase : ∀ a : ↥(unramifiedCompletionInt K),
      algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) a ∈ T)
    (hsub : torsionCompletionSet K hq hϖmax hϖne0 g hg0 hg1 hg n
      ⊆ (T : Set ↥(closureCompletionInt K))) :
    torsionCompletionSet K hq hπmax hπne0 f hf0 hf1 hf n
      ⊆ (T : Set ↥(closureCompletionInt K)) := by
  rintro w ⟨lam, hlam, hw⟩
  have hnorm := norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
    lam hlam
  obtain ⟨r, hr, hrv⟩ := exists_mem_torsionPoints_evalAt K hq hπmax hπne0 f hf0 hf1 hf
    hϖmax hϖne0 g hg0 hg1 hg σ θ hθ0 hint n lam hlam
  have hy : evalAt K θ hnorm ∈ T := hsub ⟨r, hr, hrv⟩
  have hmem : PowerSeries.aeval (hasEval_evalAt K hnorm θ hθ0) θ' ∈ T :=
    aeval_mem_of_isClosed T hT hbase (hasEval_evalAt K hnorm θ hθ0) hy θ'
  have heq : PowerSeries.aeval (hasEval_evalAt K hnorm θ hθ0) θ' = w :=
    Subtype.ext (by rw [coe_aeval_evalAt_comp_inverse K hnorm θ θ' hθ0 hθ'θ, ← hw])
  rwa [heq] at hmem

/-! ## 3. `K̂^m_f` —— 原典の主張そのもの -/

/-- **`K̂^m_f`** —— `K̂^{ur}` に `Λ_{f,n}` の(`ℂ_K` での)像を添加した中間体。

原典 `L^m_f = L(µ_{f,m})`(`L = K̂`)。★`ℂ_K` の中で取っているのは、`θ` の係数が
`𝒪_{K̂^{ur}}` にあり、`λ` が `K^al` にあるという**別々の完備化**を 1 つの体で
会わせるため(M3 の設計)。 -/
noncomputable def lubinTateCompletionField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) : IntermediateField (unramifiedCompletion K) (closureCompletion K) :=
  IntermediateField.adjoin (unramifiedCompletion K)
    (closureCompletionCoe K '' ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n))

/-- **`K̂^m_f` は `K̂^{ur}` 上有限次元**——`Λ_{f,n}` は有限集合で、各元は整
(`isIntegral_closureCompletionCoe`)だから。 -/
theorem finiteDimensional_lubinTateCompletionField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    FiniteDimensional (unramifiedCompletion K)
      ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n) := by
  haveI : Finite ↥(closureCompletionCoe K ''
      ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)) :=
    ((iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n).finite_toSet.image
      (closureCompletionCoe K)).to_subtype
  exact IntermediateField.finiteDimensional_adjoin
    (fun x hx => by obtain ⟨lam, -, rfl⟩ := hx; exact isIntegral_closureCompletionCoe K lam)

/-- ★★★★★★★★**`K̂^m_f` は `ℂ_K` の中で閉**——有限次元だから
(抽象核 `isClosed_intermediateField_of_finiteDimensional`)。

★★これが原典 Lemma 4.6 の「`[θ] ∈ 𝒪_L[[X]]` だから `[θ](µ_{f,m}) ⊂ L^m_f`」を
支えている事実である。これを落とすと `θ(λ) ∈ K̂^m_f` は言えない。 -/
theorem isClosed_lubinTateCompletionField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    IsClosed ((lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n :
      IntermediateField (unramifiedCompletion K) (closureCompletion K)) :
      Set (closureCompletion K)) := by
  haveI := finiteDimensional_lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
  exact isClosed_intermediateField_of_finiteDimensional _

/-! ### 中間体を `𝒪_{ℂ_K}` へ引き戻す

輸送核 `torsionCompletionSet_subset_of_isClosed` は `𝒪_{ℂ_K}` の部分環を相手に
するので、`ℂ_K` の中間体 `E` を `E ∩ 𝒪_{ℂ_K}` として引き戻しておく。 -/

/-- `ℂ_K` の中間体 `E` を `𝒪_{ℂ_K}` へ引き戻した部分環 `E ∩ 𝒪_{ℂ_K}`。 -/
noncomputable def intSubring (K : PAdicLocalField p)
    (E : IntermediateField (unramifiedCompletion K) (closureCompletion K)) :
    Subring ↥(closureCompletionInt K) :=
  Subring.comap (SubringClass.subtype (closureCompletionInt K)) E.toSubalgebra.toSubring

theorem mem_intSubring (K : PAdicLocalField p)
    (E : IntermediateField (unramifiedCompletion K) (closureCompletion K))
    (w : ↥(closureCompletionInt K)) :
    w ∈ intSubring K E ↔ (w : closureCompletion K) ∈ E := Iff.rfl

theorem isClosed_intSubring (K : PAdicLocalField p)
    (E : IntermediateField (unramifiedCompletion K) (closureCompletion K))
    (hE : IsClosed ((E : Set (closureCompletion K)))) :
    IsClosed ((intSubring K E : Subring ↥(closureCompletionInt K)) :
      Set ↥(closureCompletionInt K)) :=
  hE.preimage continuous_subtype_val

theorem algebraMap_mem_intSubring (K : PAdicLocalField p)
    (E : IntermediateField (unramifiedCompletion K) (closureCompletion K))
    (a : ↥(unramifiedCompletionInt K)) :
    algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) a ∈ intSubring K E := by
  rw [mem_intSubring, algebraMap_unramifiedInt_closureCompletionInt_coe]
  exact E.algebraMap_mem _

/-! ### 片側の包含、そして一致 -/

/-- ★★★★★★★★★★★★**片側の包含 `K̂^m_f ≤ K̂^m_g`**。

`K̂^m_g` は閉(`isClosed_lubinTateCompletionField`)なので、その `𝒪_{ℂ_K}` への
引き戻しも閉。あとは輸送核に投げるだけ。 -/
theorem lubinTateCompletionField_le_of_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
      ≤ lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n := by
  rw [lubinTateCompletionField, IntermediateField.adjoin_le_iff]
  rintro x ⟨lam, hlam, rfl⟩
  have hnorm := norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
    lam (Finset.mem_coe.mp hlam)
  have hsub : torsionCompletionSet K hq hϖmax hϖne0 g hg0 hg1 hg n
      ⊆ ((intSubring K (lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n) :
          Subring ↥(closureCompletionInt K)) : Set ↥(closureCompletionInt K)) := by
    rintro w ⟨r, hr, hwr⟩
    rw [SetLike.mem_coe, mem_intSubring, hwr]
    exact IntermediateField.subset_adjoin _ _ ⟨r, hr, rfl⟩
  have hmem := torsionCompletionSet_subset_of_isClosed K hq hπmax hπne0 f hf0 hf1 hf
    hϖmax hϖne0 g hg0 hg1 hg σ θ θ' hθ0 hθ'θ hint n
    (intSubring K (lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n))
    (isClosed_intSubring K _ (isClosed_lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n))
    (algebraMap_mem_intSubring K _) hsub
    (a := (⟨closureCompletionCoe K lam, (mem_closureCompletionInt K _).mpr
      (by rw [norm_closureCompletionCoe]; exact hnorm.le)⟩ : ↥(closureCompletionInt K)))
    ⟨lam, Finset.mem_coe.mp hlam, rfl⟩
  exact (mem_intSubring K _ _).mp hmem

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Corollary 4.9 前半 —— `K̂^m_f = K̂^m_{f′}`**(`θ` を仮定に置いた形)。

原典 p.9「Lemma 4.6 shows K̂^m_f = K̂^m_{f′}」。両側の包含を出している:
順方向は `θ`、逆方向は `θ′`(絡みは抽象核 `subst_intertwine_of_comp_inverse` で
裏返す)。★`θ′` の定数項は `constantCoeff_of_subst_eq_X` で自前に導く。 -/
theorem lubinTateCompletionField_eq_of_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθθ' : PowerSeries.subst θ' θ = PowerSeries.X)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
      = lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n := by
  have hθ'0 : PowerSeries.constantCoeff θ' = 0 := constantCoeff_of_subst_eq_X hθ0 hθ'θ
  have hint' : PowerSeries.subst (PowerSeries.map (baseIntHom K) g)
      (PowerSeries.map (unramGalCompletionIntHom K σ) θ')
    = PowerSeries.subst θ' (PowerSeries.map (baseIntHom K) f) :=
    subst_intertwine_of_comp_inverse _
      (constantCoeff_map_eq_zero _
        (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0))
      hθ0 hθ'0 hθθ' hθ'θ hint
  exact le_antisymm
    (lubinTateCompletionField_le_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
      hϖmax hϖne0 g hg0 hg1 hg σ θ θ' hθ0 hθ'θ hint n)
    (lubinTateCompletionField_le_of_dworkTheta K hq hϖmax hϖne0 g hg0 hg1 hg
      hπmax hπne0 f hf0 hf1 hf σ θ' θ hθ'0 hθθ' hint' n)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Corollary 4.9 前半の完成形** —— `θ` を仮定に置かず、Λ6(Dwork)から取ってくる。

`f ∈ F_π` と `g ∈ F_ϖ`(`ϖ = uπ`、`u` は単数)について `K̂^m_f = K̂^m_g`。
原典 Proposition 4.8(`θ ∈ Θ^{K̂,×}_{π,π′}` の存在)は
`exists_arithFrobenius_isCoherent_dworkThetaStep2` が供給する。

★★退化の自己検査: `u` が単数であることを落とすと `span {uπ} = span {π}` が
言えず、`Λ_{g,n}` の定義に必要な `maximalIdeal = span {uπ}` がそもそも作れない。
すなわち「同じ素元の塔」という条件は statement に埋め込まれている。 -/
theorem lubinTateCompletionField_eq (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
      = lubinTateCompletionField K hq
          (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
          (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg n := by
  obtain ⟨σ, -, -, -, hstep2⟩ := exists_arithFrobenius_isCoherent_dworkThetaStep2 K hq
  obtain ⟨θ, θ', hθ0, -, hθθ', hθ'θ, -, hint⟩ := hstep2 π hπmax f hf0 hf1 hf u hu g hg0 hg1 hg
  exact lubinTateCompletionField_eq_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
    (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
    (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg σ θ θ'
    (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hθ0) hθθ' hθ'θ hint n

/-! ## 4. 整数側の影(原典には無い。同じ輸送核から出る)

`𝒪_{K̂^{ur}}` と `Λ_{f,n}` の像が生成する `𝒪_{ℂ_K}` の**閉**部分環も `f` に依らない。
★位相閉包を取った形で述べている——`𝒪_{K̂^{ur}}[Λ_{f,n}]` 自身が閉であること
(完備離散付値環上の束の閉性)は本ファイルでは扱っていない。原典 Corollary 4.9 に
対応するのは §3(体の側)である。 -/

/-- `𝒪_{K̂^{ur}}` と `X` が生成する `𝒪_{ℂ_K}` の閉部分環。 -/
noncomputable def completionAdjoin (K : PAdicLocalField p) (X : Set ↥(closureCompletionInt K)) :
    Subring ↥(closureCompletionInt K) :=
  (Subring.closure
      (Set.range (algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K)) ∪ X)
    ).topologicalClosure

theorem isClosed_completionAdjoin (K : PAdicLocalField p) (X : Set ↥(closureCompletionInt K)) :
    IsClosed ((completionAdjoin K X : Subring ↥(closureCompletionInt K)) :
      Set ↥(closureCompletionInt K)) :=
  Subring.isClosed_topologicalClosure _

theorem algebraMap_mem_completionAdjoin (K : PAdicLocalField p)
    (X : Set ↥(closureCompletionInt K)) (a : ↥(unramifiedCompletionInt K)) :
    algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) a ∈ completionAdjoin K X :=
  Subring.le_topologicalClosure _ (Subring.subset_closure (Or.inl ⟨a, rfl⟩))

theorem subset_completionAdjoin (K : PAdicLocalField p) (X : Set ↥(closureCompletionInt K)) :
    X ⊆ (completionAdjoin K X : Set ↥(closureCompletionInt K)) :=
  fun _ hx => Subring.le_topologicalClosure _ (Subring.subset_closure (Or.inr hx))

theorem completionAdjoin_le (K : PAdicLocalField p) {X Y : Set ↥(closureCompletionInt K)}
    (h : X ⊆ (completionAdjoin K Y : Set ↥(closureCompletionInt K))) :
    completionAdjoin K X ≤ completionAdjoin K Y :=
  Subring.topologicalClosure_minimal _
    (Subring.closure_le.mpr (Set.union_subset
      (fun _ hz => by obtain ⟨a, rfl⟩ := hz; exact algebraMap_mem_completionAdjoin K Y a) h))
    (isClosed_completionAdjoin K Y)

theorem completionAdjoin_eq (K : PAdicLocalField p) {X Y : Set ↥(closureCompletionInt K)}
    (hXY : X ⊆ (completionAdjoin K Y : Set ↥(closureCompletionInt K)))
    (hYX : Y ⊆ (completionAdjoin K X : Set ↥(closureCompletionInt K))) :
    completionAdjoin K X = completionAdjoin K Y :=
  le_antisymm (completionAdjoin_le K hXY) (completionAdjoin_le K hYX)

/-- 輸送核を `T := completionAdjoin K (Λ_{g,n} の像)` に当てただけ。 -/
theorem torsionCompletionSet_subset_completionAdjoin (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    torsionCompletionSet K hq hπmax hπne0 f hf0 hf1 hf n
      ⊆ (completionAdjoin K (torsionCompletionSet K hq hϖmax hϖne0 g hg0 hg1 hg n) :
          Set ↥(closureCompletionInt K)) :=
  torsionCompletionSet_subset_of_isClosed K hq hπmax hπne0 f hf0 hf1 hf hϖmax hϖne0 g hg0 hg1 hg
    σ θ θ' hθ0 hθ'θ hint n _ (isClosed_completionAdjoin K _)
    (algebraMap_mem_completionAdjoin K _) (subset_completionAdjoin K _)

/-- `𝒪_{K̂^{ur}}[Λ_{f,n}]` の閉包は `f` に依らない(`θ` を仮定に置いた形)。 -/
theorem completionAdjoin_torsionCompletionSet_eq_of_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθθ' : PowerSeries.subst θ' θ = PowerSeries.X)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    completionAdjoin K (torsionCompletionSet K hq hπmax hπne0 f hf0 hf1 hf n)
      = completionAdjoin K (torsionCompletionSet K hq hϖmax hϖne0 g hg0 hg1 hg n) := by
  have hθ'0 : PowerSeries.constantCoeff θ' = 0 := constantCoeff_of_subst_eq_X hθ0 hθ'θ
  have hint' : PowerSeries.subst (PowerSeries.map (baseIntHom K) g)
      (PowerSeries.map (unramGalCompletionIntHom K σ) θ')
    = PowerSeries.subst θ' (PowerSeries.map (baseIntHom K) f) :=
    subst_intertwine_of_comp_inverse _
      (constantCoeff_map_eq_zero _
        (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0))
      hθ0 hθ'0 hθθ' hθ'θ hint
  exact completionAdjoin_eq K
    (torsionCompletionSet_subset_completionAdjoin K hq hπmax hπne0 f hf0 hf1 hf
      hϖmax hϖne0 g hg0 hg1 hg σ θ θ' hθ0 hθ'θ hint n)
    (torsionCompletionSet_subset_completionAdjoin K hq hϖmax hϖne0 g hg0 hg1 hg
      hπmax hπne0 f hf0 hf1 hf σ θ' θ hθ'0 hθθ' hint' n)

/-- `𝒪_{K̂^{ur}}[Λ_{f,n}]` の閉包は `f` に依らない(完成形)。 -/
theorem completionAdjoin_torsionCompletionSet_eq (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    completionAdjoin K (torsionCompletionSet K hq hπmax hπne0 f hf0 hf1 hf n)
      = completionAdjoin K (torsionCompletionSet K hq
          (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
          (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg n) := by
  obtain ⟨σ, -, -, -, hstep2⟩ := exists_arithFrobenius_isCoherent_dworkThetaStep2 K hq
  obtain ⟨θ, θ', hθ0, -, hθθ', hθ'θ, -, hint⟩ := hstep2 π hπmax f hf0 hf1 hf u hu g hg0 hg1 hg
  exact completionAdjoin_torsionCompletionSet_eq_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
    (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
    (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg σ θ θ'
    (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hθ0) hθθ' hθ'θ hint n

end ABC3.Found.PGC

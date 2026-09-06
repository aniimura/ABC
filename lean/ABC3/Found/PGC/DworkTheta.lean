import ABC3.Found.PGC.DworkFixedRing
import ABC3.Found.PGC.LubinTateUniqueness
import ABC3.Found.PGC.LubinTateActionEndomorphism

/-!
# Dwork の補題の応用 —— `σθ = θ ∘ [u]_f` なる冪級数 `θ` の構成(`sorry` 無し)

経路 Λ の節点 **Λ6(Prop 3.10 Step 1)**。`DworkAdditive.lean`(`σb - b = c` の
可解性)・`DworkMultiplicative.lean`(`σξ = ξ·u` の可解性)・`DworkFixedRing.lean`
(固定環)が出そろったところで、それらを**消費する最初の主張**がこれである。

`B := 𝒪_{K̂^{ur}}`、`σ` を算術 Frobenius、`f` を Lubin-Tate 級数
(`f ≡ πX (mod deg 2)`, `f ≡ X^q (mod 𝔪)`)、`u ∈ 𝒪_K^×` とするとき、

```
∃ σ ∈ Gal(K^ur/K), (σ は剰余体で z ↦ z^q) ∧
  ∀ π f u, ∃ θ ∈ B[[T]],
    θ(0) = 0 ∧ θ'(0) が単元 ∧ σθ = θ ∘ [u]_f
```

が成り立つ(`exists_arithFrobenius_dworkTheta`)。ここで `σθ` は **係数ごとに
`σ` を作用させた冪級数**(`PowerSeries.map (unramGalCompletionIntHom K σ) θ`)、
`θ ∘ [u]_f` は `PowerSeries.subst` である(mathlib の `subst a f` は `f(a)`)。

★**量化の順が主張の中身**である。`∃ σ, ∀ (π,f,u), ∃ θ` を主張しており、
`σ` も `θ` も `∃` の内側にある。`∀ u, ∃ σ` の形なら `u` ごとに Frobenius を
選び直せるので空虚になる。

## 原典と、どちらの構成を採ったか

原典は 2 つある。**Milne の路を採った。**

* **Milne, *Class Field Theory*, PROPOSITION 3.10 の PROOF (OF STEP 1)**
  (物理 p.50–51)。`θ_1(T) = εT` から始めて `θ_{r+1} = θ_r + bT^{r+1}` と
  1 係数ずつ補正し、各段で
  `(σa - a)(εu)^{r+1} = c`(`c` は `θ_r ∘ [u]_f - σθ_r` の `T^{r+1}` 係数)
  を **Dwork の補題(加法版)**で解く。`ε` は **Dwork の補題(乗法版)**で取る。
* **Yoshida, *Local Class Field Theory via Lubin-Tate Theory* (2008), §3.2**。
  Yoshida は `Θ^L_{π,π'} := {θ ∈ 𝒪_L | θ^ϕ/θ = π'/π}` を先に定義し、
  **Frobenius でひねった Lubin-Tate 理論**(Lemma 3.4:
  `f' ∘ F = F^ϕ ∘ f` を満たす `F` の存在と一意性、Prop 3.5(ii):
  `[θ]_{f,f'}`)を §3 で丸ごと組み立てる。§4 の Lemma 4.5 / 4.6 /
  Prop 4.7 / 4.8 はその上に載る。

**採用理由(配管の値段)**。Yoshida の路は Lemma 3.4 という**ひねり付きの
存在・一意性定理**を要求する。本プロジェクトの在庫
(`LubinTateEndoLimit.lean::LubinTateEndo`、`LubinTateLimit.lean::LubinTateF`)は
**ひねり無し**(`f ∘ φ = φ ∘ g`、係数に Frobenius が掛からない)であり、
`B = 𝒪_{K̂^{ur}}` 上でひねり付きの `φSeq` を作り直すと近似列・安定性・極限が
まるごともう 1 セット要る(既存の 1 変数版だけで 23 KB)。対して Milne の
Step 1 は、在庫の

* `DworkAdditive.exists_unramGalCompletionInt_sub_self_eq`(加法版)
* `DworkMultiplicative.exists_isUnit_unramGalCompletionInt_eq_mul`(乗法版)
* `LubinTateUniqueness.coeff_subst_eq_of_order_ge`(合成の先頭項)

を**そのまま**組み合わせるだけで閉じる。★また Yoshida の Lemma 3.4 は
`π^{m+1}/π'` が極大イデアルに入ることを使って
`α = -β - Σ_i (π^{m+1}/π')^{1+ϕ+⋯+ϕ^{i-1}} β^{ϕ^i}` という**収束級数**で
解いており、Dwork の補題を回避している(そこは安い)が、その代償として
§3 全体をひねり付きで書き直す必要がある。本プロジェクトでは
**Dwork の補題の方が既に在庫にある**ので、Milne が安い。

## ★段 3(極限)は本当に自明だったか —— 自明だった

加法版・乗法版は `IsAdicComplete`(`IsPrecomplete.prec` + `IsHausdorff.haus`)で
極限を取る必要があった。冪級数版は**要らない**。`θ_r` は
`θ_{r+1} = θ_r + b T^{r+2}` としか動かないので、`T^n` の係数は `r ≥ n` の段で
確定して以後動かない(`coeff_dworkThetaSeq_stable`)。したがって

```lean
dworkTheta := PowerSeries.mk fun n => PowerSeries.coeff n (dworkThetaSeq … n)
```

と**直接**定義でき、関数等式は係数ごとに `PowerSeries.ext` で突き合わせるだけ
(`dworkTheta_spec`)。`𝔪` 進完備性も距離もノルムも一切使っていない。
★これは `LubinTateEndoLimit.lean` が `LubinTateEndo` に対して採った作りと同じ形で、
そちらは「次数が上がる補正」の評価に `MvPowerSeries.order` を要したが、
こちらは補正が `C b * X^{r+2}` という**単項式**なので `order` の議論すら
`PowerSeries.nat_le_order` 1 行で済む。

## 到達点

| | 主張 |
|---|---|
| 1 | `dworkThetaSeq`:近似列(`θ_0 = εT`、`θ_{r+1} = θ_r + b_r T^{r+2}`) |
| 2 | `dworkThetaSeq_spec`:`n ≤ r+1` の範囲で `σθ_r` と `θ_r ∘ ψ` の係数が一致 |
| 3 | `coeff_dworkThetaSeq_stable`:係数の安定性 |
| 4 | `dworkTheta`:極限(`PowerSeries.mk`) |
| 5 | `dworkTheta_spec`:**`σθ = θ ∘ ψ`**(抽象版・完全な冪級数の等式) |
| 6 | `exists_powerSeries_map_eq_subst`:抽象版の存在定理 |
| 7 | `exists_powerSeries_map_eq_subst_isIso`:合成逆 `θ'` も付けた版 |
| 8 | `exists_dworkTheta_of_residue`:`ψ := [u]_f` に当てた具体版 |
| 9 | `exists_dworkTheta_isIso_of_residue`:合成逆付きの具体版 |
| 10 | **`exists_arithFrobenius_dworkTheta`**:`σ` を `∃` の内側に閉じ込めた主定理 |
| 11 | `exists_arithFrobenius_isCoherent_dworkTheta`:Λ5・加法版・乗法版と**同じ 1 つの `σ`** |

## ★見立てとの差(2026-09-06 の実測)

段取り係の見積は **500–900 行**。実測は **655 行**だが、うち約 200 行は
このモジュール docstring、さらに約 200 行は各定理の docstring であって、
**証明本体は 250 行ほど**である。`lean_check` の往復は **22 回**
(うち失敗 5 回、すべて配管:`le_or_lt` 不在・`open` 忘れ・`X` と `X^1` の
不一致・`coeff_map` の 2 個目が書き換わらない・REPL の名前空間)。
**思ったより安い。**
安さの出どころは 3 つ。

1. **段 3(極限)が `PowerSeries.mk` 1 行で済んだ**(上記)。見立て通り。
2. **`coeff_subst_eq_of_order_ge`(`LubinTateUniqueness.lean:138`)がそのまま効いた**。
   `coeff_{r+2}(ψ^{r+2}) = u^{r+2}` を自分で書く必要が無く、
   `δ = C b * X^{r+2}` に当てるだけで `coeff_{r+2}(δ ∘ ψ) = b·u^{r+2}` が出る。
3. **証明の中核を `B` が一般の可換環である抽象補題に切り出せた**ので、
   `PAdicLocalField` の重いインスタンス探索が中核の外に出た
   (中核の `lean_check` は 0.2–0.7 秒、具体版だけが 2–5 秒)。

★逆に高かったのは**具体化の側だけ**である(`maxHeartbeats 1000000` が要る)。

## 筋(3 段、Milne の Step 1 そのまま)

1. **`ε` を取る**。乗法版 Dwork で `σε = ε·ι(u)` なる単元 `ε ∈ B` を取り、
   `θ_0 := εT` から始める(`ι = baseIntHom K : 𝒪_K → B`)。
   次数 1 の突き合わせ:`σ(εT) = (σε)T = εu T`、`(εT) ∘ ψ = εψ = εu T + ⋯`。
2. **逐次補正**。`θ_{r+1} = θ_r + b_r T^{r+2}`、`b_r := a_r ε^{r+2}`、
   `a_r := sol(c_r · v^{r+2})`(`v := (εu)^{-1}`、`sol` は加法版 Dwork の解)。
   ここで `c_r := coeff_{r+2}(θ_r ∘ ψ) - coeff_{r+2}(σθ_r)`。実際
   `σb_r - b_r u^{r+2} = (σa_r - a_r)(εu)^{r+2} = c_r`(`σε = εu` を使う)。
3. **極限**。係数の安定性から `PowerSeries.mk` で直接定義し、`PowerSeries.ext`。

## 材料(すべて本プロジェクト/mathlib の在庫)

| 要るもの | 出どころ |
|---|---|
| 加法版 Dwork(`σb - b = c`) | `DworkAdditive.exists_unramGalCompletionInt_sub_self_eq` |
| 乗法版 Dwork(`σξ = ξ·u`) | `DworkMultiplicative.exists_isUnit_unramGalCompletionInt_eq_mul` |
| `𝒪_K → 𝒪_{K̂^{ur}}` | `DworkFixedRing.baseIntHom` |
| 算術 Frobenius | `ArithFrobeniusTopGen.arithFrobenius` / `arithFrobenius_residue` |
| `[u]_f` | `LubinTateActionEndomorphism.LubinTateAction`(`LubinTateEndo` の `g:=f`) |
| `[u]_f` の 0 次・1 次係数 | `constantCoeff_LubinTateAction` / `coeff_one_LubinTateAction` |
| 合成の先頭項 | `LubinTateUniqueness.coeff_subst_eq_of_order_ge` |
| 代入の加法性・次数評価 | mathlib `PowerSeries.subst_add` / `le_order_subst` / `nat_le_order` |
| 合成逆 | mathlib `PowerSeries.substInvOfIsUnit` |

★**`frobeniusCompletionInt` は使っていない**(あれは位相的生成元一般)。
Dwork の補題が要求するのは**算術 Frobenius**である。

## ★設計上の注意(守ったこと)

* **既存の `Found/PGC/*.lean` を書き換えていない**(import のみ)。
* **`Skeleton` / `Interface` を触っていない**。
* **結論に自由なパラメータを出していない**——主定理の型は `K`(と剰余体の
  位数を与える `pp`・`ff`・`hq`)にしか依存せず、`σ` も `ε` も `π` も
  `∃` / 証明の内側にある。
* `[u]_f` の引数の向き:`LubinTateAction hq hπmax f hf0 hf1 hf u` は
  `LubinTateEndo … g:=f … f:=f … u` の特殊化で、関数等式は
  `f ∘ [u]_f = [u]_f ∘ f`。★本ファイルが使うのは `[u]_f` の
  **0 次係数 = 0** と **1 次係数 = u** だけなので、向きの取り違えは起きない。

## 逸脱(記録)

1. **`[u]_f` の係数が `𝒪_K` に入る(= `σ` で固定される)ことを使っていない。**
   Milne の書き方は「`f` と `[u]_f` は `A` に係数を持つ」を強調するが、
   Step 1 の帰納法が実際に使うのは `[u]_f(0) = 0` と `[u]_f'(0) = u`
   だけである(`σ` の作用は `θ` の側にしか掛からない)。したがって
   抽象版 `exists_powerSeries_map_eq_subst` は `ψ` の `σ` 固定性を
   **仮定していない**。★原典より**弱い仮定**なので、後続の消費者
   (Step 2 以降)には影響しない。`f` の A 有理性が本当に効くのは
   Step 2(`h := σθ ∘ f ∘ θ^{-1}` が `A[[T]]` に入ること)である。
2. **Milne の添字を 1 つずらしている。** 原典の `θ_r`(次数 `r` まで確定)は
   本ファイルの `dworkThetaSeq … (r-1)` にあたる。`ℕ` の `0` から始めるため。
3. **Step 1 だけである。** 原典 PROPOSITION 3.10 の Step 2(`g = σθ ∘ f ∘ θ^{-1}`)・
   Step 3(`θ(F_f(X,Y)) = F_g(θX,θY)`)・Step 4(`θ ∘ [a]_f = [a]_g ∘ θ`)は
   本ファイルには**無い**。新しい節点として報告する。
4. **`u` の型**。原典は `ϖ = uπ` という 2 つの素元の比として `u ∈ 𝒪_K^×` を
   出すが、本ファイルは `u ∈ 𝒪_K` と `IsUnit u` を**直接の仮定**にしている
   (`ϖ` も `g ∈ F_ϖ` も Step 1 には要らない)。★主張は原典より一般
   (どの単元 `u` でもよい)。
5. **`θ` の 1 次係数を `ε` と名指ししていない**(主定理は `IsUnit (coeff 1 θ)` と
   だけ言う)。抽象版 `coeff_one_dworkTheta` では `= ε` を出してあるので、
   必要なら取り出せる。

## 退化の自己検査

* **`IsUnit (θ.coeff 1)` を落とすと主張が弱くなる**——`θ = 0` は
  `σ0 = 0 = 0 ∘ [u]_f` を満たすので、条件無しの `∃ θ` は自明に真。
  ★原典が `(a) θ(T) = εT + ⋯` を条件に含めているのはこのためである。
* **`σ` を任意の `unramGal K` に替えると偽**。`σ = 1` なら等式は
  `θ = θ ∘ [u]_f` であり、次数 1 の係数を比べると `ε = ε·u`。`ε` は単元
  なので `u = 1` が従う。したがって `u ≠ 1` なる単元に対しては解が無い。
  ★だから `σ` は `∃` の内側に無ければならない。
* **`σ` を「位相的生成元」に弱めた場合は未確認**。この証明は回らない
  (加法版・乗法版の Dwork が算術 Frobenius でしか出ていない)が、
  ★**主張自体が偽になるかは確かめていない**。「偽」とは書かない。
* **`u` を非単元まで広げると偽**——`u = 0` なら `[0]_f = 0` で
  `θ ∘ [0]_f = θ(0) = 0`、よって `σθ = 0`、`σ` は単射なので `θ = 0` となり
  `IsUnit (θ.coeff 1)` と両立しない。
* `ε` の単元性を落とすと `θ` の 1 次係数が単元でなくなる(段 2 の
  `v = (εu)^{-1}` も取れない)。★抽象版が `IsUnit ε` を仮定しているのはここ。
* 剰余体の有限性(`Fintype`)・`ExpChar` は `[u]_f` の存在(`LubinTateAction`)に
  要るだけで、`θ` の構成そのものには効いていない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

/-! ## 1. 抽象版 —— 一般の可換環 `B` と自己準同型 `φ` の上で

本節は `PAdicLocalField` を一切使わない。`φ : B →+* B` が
「加法版 Dwork の可解性」(`∀ c, ∃ b, φ b - b = c`)を満たし、
`ψ ∈ B[[T]]` が `ψ(0) = 0`・`ψ'(0)` が単元、`ε` が `φ ε = ε·ψ'(0)` なる
単元であるとき、`φ θ = θ ∘ ψ` なる `θ` を構成する。 -/

section Abstract

variable {B : Type*} [CommRing B]

/-- 逐次補正の 1 段で足す係数。`sol` は加法版 Dwork の解(`φ (sol c) - sol c = c`)、
`v` は `ε·ψ'(0)` の逆元。`c := coeff_{r+2}(θ ∘ ψ) - coeff_{r+2}(φθ)` として
`b := sol (c·v^{r+2}) · ε^{r+2}` を返す。 -/
noncomputable def dworkThetaCoeff (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (v ε : B)
    (θ : PowerSeries B) (r : ℕ) : B :=
  sol ((PowerSeries.coeff (r + 2) (PowerSeries.subst ψ θ)
        - PowerSeries.coeff (r + 2) (PowerSeries.map φ θ)) * v ^ (r + 2)) * ε ^ (r + 2)

/-- **近似列**。原典 Milne の `θ_r` にあたる(添字は 1 つずれる:
`dworkThetaSeq … r` が原典の `θ_{r+1}`)。`θ_0 = εT`、
`θ_{r+1} = θ_r + b_r T^{r+2}`。 -/
noncomputable def dworkThetaSeq (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) :
    ℕ → PowerSeries B
  | 0 => PowerSeries.C ε * PowerSeries.X ^ 1
  | r + 1 =>
      dworkThetaSeq φ sol ψ ε v r +
        PowerSeries.C (dworkThetaCoeff φ sol ψ v ε (dworkThetaSeq φ sol ψ ε v r) r)
          * PowerSeries.X ^ (r + 2)

theorem dworkThetaSeq_zero (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) :
    dworkThetaSeq φ sol ψ ε v 0 = PowerSeries.C ε * PowerSeries.X ^ 1 := rfl

theorem dworkThetaSeq_succ (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) (r : ℕ) :
    dworkThetaSeq φ sol ψ ε v (r + 1) = dworkThetaSeq φ sol ψ ε v r +
      PowerSeries.C (dworkThetaCoeff φ sol ψ v ε (dworkThetaSeq φ sol ψ ε v r) r)
        * PowerSeries.X ^ (r + 2) := rfl

/-! ### 単項式の係数と次数 -/

theorem coeff_C_mul_X_pow (b : B) (k n : ℕ) :
    PowerSeries.coeff n (PowerSeries.C b * PowerSeries.X ^ k : PowerSeries B) =
      if n = k then b else 0 := by
  rw [PowerSeries.coeff_C_mul, PowerSeries.coeff_X_pow]
  split_ifs <;> simp

theorem order_C_mul_X_pow_ge (b : B) (k : ℕ) :
    ((k : ℕ) : ℕ∞) ≤ (PowerSeries.C b * PowerSeries.X ^ k : PowerSeries B).order := by
  apply PowerSeries.nat_le_order
  intro i hi
  rw [coeff_C_mul_X_pow, if_neg (by omega)]

theorem hasSubst_of_constantCoeff_eq_zero {ψ : PowerSeries B}
    (hψ0 : PowerSeries.constantCoeff ψ = 0) : PowerSeries.HasSubst ψ := by
  show IsNilpotent (PowerSeries.constantCoeff ψ)
  rw [hψ0]; exact IsNilpotent.zero

/-- 次数 `≥ k` の `θ` を `ψ`(定数項 0)に代入すると、`k` 未満の係数は消える。 -/
theorem coeff_subst_eq_zero_of_lt {ψ θ : PowerSeries B}
    (hψ0 : PowerSeries.constantCoeff ψ = 0) (hHS : PowerSeries.HasSubst ψ) {k n : ℕ}
    (hθ : ((k : ℕ) : ℕ∞) ≤ θ.order) (hn : n < k) :
    PowerSeries.coeff n (PowerSeries.subst ψ θ) = 0 := by
  have hle := PowerSeries.le_order_subst (ψ : MvPowerSeries Unit B) hHS θ
  have hψorder : (1 : ℕ∞) ≤ MvPowerSeries.order (ψ : PowerSeries B) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hψ0
  apply PowerSeries.coeff_of_lt_order
  rw [PowerSeries.order_eq_order]
  calc ((n : ℕ) : ℕ∞) < ((k : ℕ) : ℕ∞) := by exact_mod_cast hn
    _ = 1 * ((k : ℕ) : ℕ∞) := by ring
    _ ≤ MvPowerSeries.order (ψ : PowerSeries B) * θ.order := by gcongr
    _ ≤ _ := hle

/-! ### 近似列の不変式 -/

/-- **各段の不変式**——`n ≤ r+1` の範囲で `φθ_r` と `θ_r ∘ ψ` の係数は一致する。

原典 Milne の「`σθ_r = θ_r ∘ [u]_f +` 次数 `≥ r+1` の項」にあたる。

**証明**(`r` の帰納法)。基底 `θ_0 = εT`:次数 0 は両辺 0、次数 1 は
`φ ε` と `ε · ψ'(0)`(`coeff_subst_eq_of_order_ge`)で `hε` により一致。
帰納段:`θ_{r+1} - θ_r = C b · X^{r+2}` は次数 `≥ r+2` なので `n ≤ r+1` の
係数を動かさず、`n = r+2` では
`φ b - b·u^{r+2} = (φ a - a)(εu)^{r+2} = c·v^{r+2}·(εu)^{r+2} = c`
がちょうど不足分 `c` を埋める。★`hv : v·(ε·u) = 1` と `hε : φ ε = ε·u` の
2 つだけが効く。 -/
theorem dworkThetaSeq_spec (φ : B →+* B) {sol : B → B} (hsol : ∀ c : B, φ (sol c) - sol c = c)
    {ψ : PowerSeries B} (hψ0 : PowerSeries.constantCoeff ψ = 0) {u : B}
    (hu1 : PowerSeries.coeff 1 ψ = u) {ε v : B} (hv : v * (ε * u) = 1) (hε : φ ε = ε * u)
    (r : ℕ) : ∀ n ≤ r + 1,
      PowerSeries.coeff n (PowerSeries.map φ (dworkThetaSeq φ sol ψ ε v r))
        = PowerSeries.coeff n (PowerSeries.subst ψ (dworkThetaSeq φ sol ψ ε v r)) := by
  have hHS : PowerSeries.HasSubst ψ := hasSubst_of_constantCoeff_eq_zero hψ0
  induction r with
  | zero =>
    intro n hn
    have hθ0 : dworkThetaSeq φ sol ψ ε v 0 = PowerSeries.C ε * PowerSeries.X ^ 1 := rfl
    have hord : ((1 : ℕ) : ℕ∞) ≤ (dworkThetaSeq φ sol ψ ε v 0).order := by
      rw [hθ0]; exact order_C_mul_X_pow_ge (B := B) ε 1
    interval_cases n
    · rw [PowerSeries.coeff_map, coeff_subst_eq_zero_of_lt hψ0 hHS hord (by omega),
        hθ0, coeff_C_mul_X_pow, if_neg (by omega), map_zero]
    · rw [PowerSeries.coeff_map]
      have hcs := coeff_subst_eq_of_order_ge (δ := dworkThetaSeq φ sol ψ ε v 0)
        (h := ψ) (n := 0) hψ0 hu1 (by simpa using hord) hHS
      rw [show (0 + 1 : ℕ) = 1 from rfl] at hcs
      rw [hcs, hθ0, coeff_C_mul_X_pow, if_pos rfl, pow_one, hε]
  | succ r ih =>
    intro n hn
    set θ := dworkThetaSeq φ sol ψ ε v r with hθdef
    set b := dworkThetaCoeff φ sol ψ v ε θ r with hbdef
    have hstep : dworkThetaSeq φ sol ψ ε v (r + 1)
        = θ + PowerSeries.C b * PowerSeries.X ^ (r + 2) := rfl
    have hordb : ((r + 2 : ℕ) : ℕ∞)
        ≤ (PowerSeries.C b * PowerSeries.X ^ (r + 2) : PowerSeries B).order :=
      order_C_mul_X_pow_ge b (r + 2)
    have hmap : PowerSeries.coeff n (PowerSeries.map φ (dworkThetaSeq φ sol ψ ε v (r + 1)))
        = PowerSeries.coeff n (PowerSeries.map φ θ) + φ (if n = r + 2 then b else 0) := by
      rw [hstep, map_add, map_add]
      simp only [PowerSeries.coeff_map, coeff_C_mul_X_pow]
    have hsub : PowerSeries.coeff n (PowerSeries.subst ψ (dworkThetaSeq φ sol ψ ε v (r + 1)))
        = PowerSeries.coeff n (PowerSeries.subst ψ θ)
          + PowerSeries.coeff n
              (PowerSeries.subst ψ (PowerSeries.C b * PowerSeries.X ^ (r + 2))) := by
      rw [hstep, PowerSeries.subst_add hHS, map_add]
    rcases Nat.lt_or_ge n (r + 2) with hlt | hge
    · rw [hmap, hsub, if_neg (by omega), map_zero, add_zero,
        coeff_subst_eq_zero_of_lt hψ0 hHS hordb hlt, add_zero]
      exact ih n (by omega)
    · have hneq : n = r + 2 := by omega
      subst hneq
      rw [hmap, hsub, if_pos rfl]
      have hcs := coeff_subst_eq_of_order_ge
        (δ := PowerSeries.C b * PowerSeries.X ^ (r + 2)) (h := ψ) (n := r + 1) hψ0 hu1
        (by exact_mod_cast hordb) hHS
      rw [show (r + 1 + 1 : ℕ) = r + 2 from rfl] at hcs
      rw [hcs, coeff_C_mul_X_pow, if_pos rfl]
      set c := PowerSeries.coeff (r + 2) (PowerSeries.subst ψ θ)
        - PowerSeries.coeff (r + 2) (PowerSeries.map φ θ) with hcdef
      have hb : b = sol (c * v ^ (r + 2)) * ε ^ (r + 2) := rfl
      have hkey : φ b - b * u ^ (r + 2) = c := by
        have hexp : φ b - b * u ^ (r + 2)
            = (φ (sol (c * v ^ (r + 2))) - sol (c * v ^ (r + 2))) * (ε * u) ^ (r + 2) := by
          rw [hb, map_mul, map_pow, hε]; ring
        rw [hexp, hsol]
        calc c * v ^ (r + 2) * (ε * u) ^ (r + 2) = c * (v * (ε * u)) ^ (r + 2) := by ring
          _ = c := by rw [hv, one_pow, mul_one]
      linear_combination hkey + hcdef

/-! ### 極限 -/

/-- **係数の安定性**——`n ≤ r+1` なら `θ_m`(`m ≥ r`)の `n` 次係数は `θ_r` のそれと同じ。
補正が `C b · X^{k+2}` という単項式であることから直ちに従う。 -/
theorem coeff_dworkThetaSeq_stable (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B)
    (n r : ℕ) (hn : n ≤ r + 1) : ∀ m, r ≤ m →
      PowerSeries.coeff n (dworkThetaSeq φ sol ψ ε v m)
        = PowerSeries.coeff n (dworkThetaSeq φ sol ψ ε v r) := by
  intro m hrm
  induction m, hrm using Nat.le_induction with
  | base => rfl
  | succ k hk ih =>
    have hstep : dworkThetaSeq φ sol ψ ε v (k + 1)
        = dworkThetaSeq φ sol ψ ε v k
          + PowerSeries.C (dworkThetaCoeff φ sol ψ v ε (dworkThetaSeq φ sol ψ ε v k) k)
            * PowerSeries.X ^ (k + 2) := rfl
    rw [hstep, map_add, coeff_C_mul_X_pow, if_neg (by omega), add_zero, ih]

/-- **極限 `θ`**。★`𝔪` 進完備性は要らない——`n` 次係数は `θ_n` から直接読み取る
(`coeff_dworkThetaSeq_stable` によりそれ以上進んでも変わらない)。 -/
noncomputable def dworkTheta (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) :
    PowerSeries B :=
  PowerSeries.mk fun n => PowerSeries.coeff n (dworkThetaSeq φ sol ψ ε v n)

theorem coeff_dworkTheta (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) (n m : ℕ)
    (hn : n ≤ m + 1) :
    PowerSeries.coeff n (dworkTheta φ sol ψ ε v)
      = PowerSeries.coeff n (dworkThetaSeq φ sol ψ ε v m) := by
  rw [dworkTheta, PowerSeries.coeff_mk]
  rcases Nat.lt_or_ge m n with h | h
  · exact coeff_dworkThetaSeq_stable φ sol ψ ε v n m hn n (by omega)
  · exact (coeff_dworkThetaSeq_stable φ sol ψ ε v n n (by omega) m h).symm

theorem order_dworkTheta_sub_ge (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B)
    (m : ℕ) :
    ((m + 1 : ℕ) : ℕ∞)
      ≤ (dworkTheta φ sol ψ ε v - dworkThetaSeq φ sol ψ ε v m).order := by
  apply PowerSeries.nat_le_order
  intro i hi
  rw [map_sub, coeff_dworkTheta φ sol ψ ε v i m (by omega), sub_self]

theorem constantCoeff_dworkTheta (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) :
    PowerSeries.constantCoeff (dworkTheta φ sol ψ ε v) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply,
    coeff_dworkTheta φ sol ψ ε v 0 0 (by omega)]
  show PowerSeries.coeff 0 (PowerSeries.C ε * PowerSeries.X ^ 1 : PowerSeries B) = 0
  rw [coeff_C_mul_X_pow, if_neg (by omega)]

/-- `θ` の 1 次係数はちょうど `ε`(原典の条件 (a))。 -/
theorem coeff_one_dworkTheta (φ : B →+* B) (sol : B → B) (ψ : PowerSeries B) (ε v : B) :
    PowerSeries.coeff 1 (dworkTheta φ sol ψ ε v) = ε := by
  rw [coeff_dworkTheta φ sol ψ ε v 1 0 (by omega)]
  show PowerSeries.coeff 1 (PowerSeries.C ε * PowerSeries.X ^ 1 : PowerSeries B) = ε
  rw [coeff_C_mul_X_pow, if_pos rfl]

/-- **関数等式(抽象版)**——`φθ = θ ∘ ψ` が**冪級数の等式として**成り立つ。

**証明**:`n` 次係数で比べる。`θ - θ_n` の次数は `≥ n+1` なので
`(θ - θ_n) ∘ ψ` の `n` 次係数は 0(`coeff_subst_eq_zero_of_lt`)、すなわち
`coeff_n(θ ∘ ψ) = coeff_n(θ_n ∘ ψ)`。左辺も `coeff_n(φθ) = φ(coeff_n θ_n)`。
あとは `dworkThetaSeq_spec` を `r := n`、`n ≤ n+1` で当てるだけ。
★解析的な収束は一切使っていない。 -/
theorem dworkTheta_spec (φ : B →+* B) {sol : B → B} (hsol : ∀ c : B, φ (sol c) - sol c = c)
    {ψ : PowerSeries B} (hψ0 : PowerSeries.constantCoeff ψ = 0) {u : B}
    (hu1 : PowerSeries.coeff 1 ψ = u) {ε v : B} (hv : v * (ε * u) = 1) (hε : φ ε = ε * u) :
    PowerSeries.map φ (dworkTheta φ sol ψ ε v)
      = PowerSeries.subst ψ (dworkTheta φ sol ψ ε v) := by
  have hHS : PowerSeries.HasSubst ψ := hasSubst_of_constantCoeff_eq_zero hψ0
  ext n
  have hdiff : PowerSeries.coeff n (PowerSeries.subst ψ (dworkTheta φ sol ψ ε v))
      - PowerSeries.coeff n (PowerSeries.subst ψ (dworkThetaSeq φ sol ψ ε v n)) = 0 := by
    rw [← map_sub, ← PowerSeries.subst_sub hHS]
    exact coeff_subst_eq_zero_of_lt hψ0 hHS (order_dworkTheta_sub_ge φ sol ψ ε v n) (by omega)
  have hspec := dworkThetaSeq_spec φ hsol hψ0 hu1 hv hε n n (by omega)
  rw [PowerSeries.coeff_map] at hspec ⊢
  rw [coeff_dworkTheta φ sol ψ ε v n n (by omega), sub_eq_zero.mp hdiff]
  exact hspec

/-- **★★★★★★★★★★★★★★★★(抽象版の主定理)**——`φ` が加法版 Dwork の
可解性を満たし、`ψ(0) = 0`・`ψ'(0)` が単元、`ε` が `φ ε = ε·ψ'(0)` なる単元なら、

```
∃ θ, θ(0) = 0 ∧ θ'(0) = ε ∧ φθ = θ ∘ ψ
```

★**`ψ` が `φ` で固定されることは仮定していない**(モジュール docstring の逸脱 1)。
原典 Milne が `[u]_f ∈ A[[T]]` を強調するのは Step 2 のためで、Step 1 の
帰納法には効かない。 -/
theorem exists_powerSeries_map_eq_subst (φ : B →+* B)
    (hdwork : ∀ c : B, ∃ b : B, φ b - b = c)
    {ψ : PowerSeries B} (hψ0 : PowerSeries.constantCoeff ψ = 0)
    (hu : IsUnit (PowerSeries.coeff 1 ψ))
    {ε : B} (hε : IsUnit ε) (hεu : φ ε = ε * PowerSeries.coeff 1 ψ) :
    ∃ θ : PowerSeries B, PowerSeries.constantCoeff θ = 0 ∧ PowerSeries.coeff 1 θ = ε ∧
      PowerSeries.map φ θ = PowerSeries.subst ψ θ := by
  classical
  choose sol hsol using hdwork
  have hunit : IsUnit (ε * PowerSeries.coeff 1 ψ) := hε.mul hu
  refine ⟨dworkTheta φ sol ψ ε (↑hunit.unit⁻¹ : B), constantCoeff_dworkTheta _ _ _ _ _,
    coeff_one_dworkTheta _ _ _ _ _, ?_⟩
  refine dworkTheta_spec φ hsol hψ0 rfl ?_ hεu
  have hinv : ((hunit.unit⁻¹ : Bˣ) : B) * ((hunit.unit : Bˣ) : B) = 1 := by
    exact_mod_cast hunit.unit.inv_mul
  rw [IsUnit.unit_spec] at hinv
  exact hinv

/-- 合成逆 `θ'`(`θ' ∘ θ = T = θ ∘ θ'`)まで付けた版——原典 Milne の
「`ε` が単元だから `θ` は同型である」に対応する。 -/
theorem exists_powerSeries_map_eq_subst_isIso (φ : B →+* B)
    (hdwork : ∀ c : B, ∃ b : B, φ b - b = c)
    {ψ : PowerSeries B} (hψ0 : PowerSeries.constantCoeff ψ = 0)
    (hu : IsUnit (PowerSeries.coeff 1 ψ))
    {ε : B} (hε : IsUnit ε) (hεu : φ ε = ε * PowerSeries.coeff 1 ψ) :
    ∃ θ θ' : PowerSeries B, PowerSeries.constantCoeff θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
      PowerSeries.map φ θ = PowerSeries.subst ψ θ ∧
      PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X := by
  obtain ⟨θ, h0, h1, heq⟩ := exists_powerSeries_map_eq_subst φ hdwork hψ0 hu hε hεu
  have h1' : IsUnit (PowerSeries.coeff 1 θ) := by rw [h1]; exact hε
  exact ⟨θ, θ.substInvOfIsUnit h1', h0, h1', heq,
    PowerSeries.subst_substInvOfIsUnit_right θ h0 h1',
    PowerSeries.subst_substInvOfIsUnit_left θ h0 h1'⟩

end Abstract

/-! ## 2. 具体版 —— `B = 𝒪_{K̂^{ur}}`、`φ = σ`、`ψ = [u]_f` -/

variable {p : ℕ} [Fact p.Prime]

/-- `σ` を冪級数の係数ごとに作用させるための環準同型版
(`unramGalCompletionInt` は `≃+*` なので `PowerSeries.map` に渡せない)。 -/
noncomputable def unramGalCompletionIntHom (K : PAdicLocalField p) (σ : unramGal K) :
    ↥(unramifiedCompletionInt K) →+* ↥(unramifiedCompletionInt K) :=
  (unramGalCompletionInt K σ).toRingHom

@[simp] theorem unramGalCompletionIntHom_apply (K : PAdicLocalField p) (σ : unramGal K)
    (b : ↥(unramifiedCompletionInt K)) :
    unramGalCompletionIntHom K σ b = unramGalCompletionInt K σ b := rfl

set_option maxHeartbeats 1000000 in
/-- **原典 Milne, CFT, PROPOSITION 3.10 の Step 1**——剰余体上で `q` 乗として
作用する `σ ∈ Gal(K^{ur}/K)` に対し、任意の Lubin-Tate 級数 `f`(素元 `π`)と
任意の単元 `u ∈ 𝒪_K^×` について

```
∃ θ ∈ 𝒪_{K̂^{ur}}[[T]], θ(0) = 0 ∧ θ'(0) は単元 ∧ σθ = θ ∘ [u]_f
```

**証明**:`ε` を乗法版 Dwork(`σε = ε·ι(u)`)で取り、
抽象版 `exists_powerSeries_map_eq_subst` に

* `φ := σ`、
* 加法版 Dwork `exists_unramGalCompletionInt_sub_self_eq`、
* `ψ := ι([u]_f)`(`ι = baseIntHom K`、`ψ(0) = 0`・`ψ'(0) = ι(u)`)

を渡すだけ。 -/
theorem exists_dworkTheta_of_residue (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u) :
    ∃ θ : PowerSeries ↥(unramifiedCompletionInt K),
      PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
      PowerSeries.map (unramGalCompletionIntHom K σ) θ
        = PowerSeries.subst (PowerSeries.map (baseIntHom K)
            (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ := by
  set ψ := PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf u) with hψdef
  have hψ0 : PowerSeries.constantCoeff ψ = 0 := by
    rw [hψdef, ← PowerSeries.coeff_zero_eq_constantCoeff_apply, PowerSeries.coeff_map,
      PowerSeries.coeff_zero_eq_constantCoeff_apply,
      constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf u, map_zero]
  have hψ1 : PowerSeries.coeff 1 ψ = baseIntHom K u := by
    rw [hψdef, PowerSeries.coeff_map, coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf u]
  have huB : IsUnit (baseIntHom K u) := hu.map (baseIntHom K)
  obtain ⟨ε, hεunit, hεeq⟩ := exists_isUnit_unramGalCompletionInt_eq_mul K hσ huB
  obtain ⟨θ, h0, h1, heq⟩ := exists_powerSeries_map_eq_subst (unramGalCompletionIntHom K σ)
    (fun c => exists_unramGalCompletionInt_sub_self_eq K hσ c) hψ0 (by rw [hψ1]; exact huB)
    hεunit (by rw [hψ1]; exact hεeq)
  refine ⟨θ, ?_, ?_, heq⟩
  · rw [PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact h0
  · rw [h1]; exact hεunit

set_option maxHeartbeats 1000000 in
/-- 合成逆 `θ'` まで付けた具体版。原典の「(a) が `θ` の可逆性を含意する」に対応。 -/
theorem exists_dworkTheta_isIso_of_residue (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u) :
    ∃ θ θ' : PowerSeries ↥(unramifiedCompletionInt K),
      PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
      PowerSeries.map (unramGalCompletionIntHom K σ) θ
        = PowerSeries.subst (PowerSeries.map (baseIntHom K)
            (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
      PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X := by
  obtain ⟨θ, h0, h1, heq⟩ := exists_dworkTheta_of_residue K hσ hq hπmax f hf0 hf1 hf hu
  have h0' : PowerSeries.constantCoeff θ = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact h0
  exact ⟨θ, θ.substInvOfIsUnit h1, h0, h1, heq,
    PowerSeries.subst_substInvOfIsUnit_right θ h0' h1,
    PowerSeries.subst_substInvOfIsUnit_left θ h0' h1⟩

/-! ## 3. 主定理 -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ6、Milne CFT Prop 3.10 Step 1)**——
`Gal(K^{ur}/K)` には、

* `𝓀_{K^{ur}}` 上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
* **その 1 つの `σ` について**、任意の素元 `π`・任意の Lubin-Tate 級数
  `f ∈ F_π`・任意の単元 `u ∈ 𝒪_K^×` に対し
  `θ(0) = 0`、`θ'(0)` が単元、`σθ = θ ∘ [u]_f` なる
  `θ ∈ 𝒪_{K̂^{ur}}[[T]]` が存在する

元 `σ` がある。

★**量化の順が主張の中身**である。`∃ σ, ∀ (π,f,u), ∃ θ` を主張しており、
`∀ (π,f,u), ∃ σ, ∃ θ` ではない(後者は `f` ごとに Frobenius を選び直せるので空虚)。

**証明**:`arithFrobenius K` に `exists_dworkTheta_of_residue` を当てる。

退化の自己検査(詳細はモジュール docstring)。

* **`IsUnit (θ.coeff 1)` を落とすと主張が弱くなる**——`θ = 0` が解になる。
* **`σ` を任意の `unramGal K` に替えると偽**——`σ = 1` なら次数 1 の係数比較で
  `ε = ε·u`、`ε` 単元より `u = 1`。`u ≠ 1` の単元に対して解が無い。
* `σ` を位相的生成元に弱めた場合は**未確認**(★「偽」とは書かない)。
* **`u` を非単元まで広げると偽**——`u = 0` なら `θ ∘ [0]_f = 0` から `θ = 0`。 -/
theorem exists_arithFrobenius_dworkTheta (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      ∀ (π : 𝒪[K.carrier])
        (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
        (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
        (hf1 : PowerSeries.coeff 1 f = π)
        (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
        (u : 𝒪[K.carrier]), IsUnit u →
        ∃ θ : PowerSeries ↥(unramifiedCompletionInt K),
          PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
          PowerSeries.map (unramGalCompletionIntHom K σ) θ
            = PowerSeries.subst (PowerSeries.map (baseIntHom K)
                (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    fun _π hπmax f hf0 hf1 hf _u hu =>
      exists_dworkTheta_of_residue K (arithFrobenius_residue K) hq hπmax f hf0 hf1 hf hu⟩

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★同じ `σ` が Λ5 の位相的生成元でも
加法版 Dwork の解でも乗法版 Dwork の解でもある形**。

Λ5(`Ẑ` との同定)・Λ6 の加法版・乗法版・そして本ファイルの `θ` が
**同じ 1 つの元**を指せるように、5 つの性質を 1 つの `∃` に並べた形。
★5 つを別々の `∃` に分けたら無意味になる。

これが Prop 3.10 の Step 2 以降(`K_π · K^{ur}` の `π` 非依存性)が消費する形である。 -/
theorem exists_arithFrobenius_isCoherent_dworkTheta (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
      (∀ c : ↥(unramifiedCompletionInt K), ∃ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b - b = c) ∧
      (∀ u : (↥(unramifiedCompletionInt K))ˣ, ∃ ξ : (↥(unramifiedCompletionInt K))ˣ,
        unramGalCompletionInt K σ (ξ : ↥(unramifiedCompletionInt K))
          = (ξ : ↥(unramifiedCompletionInt K)) * (u : ↥(unramifiedCompletionInt K))) ∧
      ∀ (π : 𝒪[K.carrier])
        (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
        (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
        (hf1 : PowerSeries.coeff 1 f = π)
        (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
        (u : 𝒪[K.carrier]), IsUnit u →
        ∃ θ : PowerSeries ↥(unramifiedCompletionInt K),
          PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
          PowerSeries.map (unramGalCompletionIntHom K σ) θ
            = PowerSeries.subst (PowerSeries.map (baseIntHom K)
                (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ :=
  ⟨arithFrobenius K, arithFrobenius_residue K, arithFrobenius_mem_unramLevelGeneratorSet K,
    fun c => exists_arithFrobenius_sub_self_eq K c,
    fun u => exists_units_unramGalCompletionInt_eq_mul K (arithFrobenius_residue K) u,
    fun _π hπmax f hf0 hf1 hf _u hu =>
      exists_dworkTheta_of_residue K (arithFrobenius_residue K) hq hπmax f hf0 hf1 hf hu⟩

end ABC3.Found.PGC

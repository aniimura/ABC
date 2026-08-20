import ABC3.Skeleton.GaloisRep.WeilRoot
import ABC3.Found.GaloisRep.WeilCharZero

/-!
# スケルトン —— **Weil 対 `e_n` とその 5 性質**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★対そのものと 3 性質が **`Found` から来た**(2026-08-20)

| 節点 | 状態 |
|---|---|
| `weilPairing`(定義) | ✅ **第 178 `weilPairingVal`** |
| `weilPairing_pow_eq_one` | ✅ **第 179** |
| `weilPairing_add_left` | ✅ **第 184 + 第 191**(`O` と `P₂ = −P₁` の場合込み) |
| `weilPairing_self` | ✅ **第 190 + 第 191** |
| `weilPairing_nondeg` | ❌ 未 |
| `weilPairing_galois` | ❌ 未 |

★対は `e_n(P, Q) := τ_Q(g_P) / g_P` で定める(`WeilRoot.lean` の `g_P`)。
★★値が `F` に落ちるのは `n` 乗が 1 だからで、第 176-178 で構成した。

### ★逸脱(記録)

3 つの仮定を足した。★どれも最終消費者 `det_galRep_eq_cyclotomic`
(`[CharZero K]`・`[IsAlgClosed L]`)が満たすので後続に影響しない。

| 仮定 | 理由 |
|---|---|
| `[W.IsElliptic]` | `weilPairingVal` の構成(群法則・生成点)が要求する |
| `[IsAlgClosed F]` | 第 137(Dedekind 性)・第 139(定数の `n` 乗根) |
| `[CharZero F]` | `μ` の存在(第 125)・`E[n²]` の位数(第 186)・`IsUnit (2:F)` |

## ★★★残る 2 性質

★★★このうち **Galois 同変性が `det ρ = 円分指標` を出す**——
`σ(e_n(P,Q)) = e_n(σP, σQ) = e_n(P,Q)^{det ρ(σ)}` が `∧²T_l E ≅ ℤ_l(1)` の内容である。
★★非退化性は `F(E)/[n]^*F(E)` の Galois 理論(Kummer 型)が要る。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

/-! ## ★★★★★対そのもの -/

/-- ★★★★★**Weil 対** `e_n : E[n] × E[n] → F`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`WeilRoot.lean` の `g_P` から `e_n(P,Q) = τ_Q(g_P)/g_P` で定める。
★★★**第 178 ブロック `weilPairingVal` がその構成である**(2026-08-20)。 -/
noncomputable def weilPairing (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) :
    W.Point → W.Point → F :=
  ABC3.Found.GaloisRep.weilPairingVal W n

/-! ## ★★5 性質 -/

/-- ★★**`e_n(P,Q)` は 1 の `n` 乗根である**。★★★第 179 ブロックで証明された。 -/
theorem weilPairing_pow_eq_one [IsAlgClosed F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (n : ℕ) (hn : 1 ≤ n) (P Q : W.Point) (hP : n • P = 0) (hQ : n • Q = 0) :
    (weilPairing W n P Q) ^ n = 1 :=
  ABC3.Found.GaloisRep.weilPairingVal_pow_eq_one W P Q hQ

/-- ★★★**双線型性**(第 1 変数)。★★★第 184 + 第 191 ブロックで証明された。 -/
theorem weilPairing_add_left [CharZero F] [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (P₁ P₂ Q : W.Point) (hP₁ : n • P₁ = 0) (hP₂ : n • P₂ = 0) (hQ : n • Q = 0) :
    weilPairing W n (P₁ + P₂) Q = weilPairing W n P₁ Q * weilPairing W n P₂ Q :=
  ABC3.Found.GaloisRep.weilPairingVal_add_left W n hn P₁ P₂ Q hP₁ hP₂ hQ

/-- ★★**交代性** `e_n(P,P) = 1`。★★★第 190 + 第 191 ブロックで証明された。 -/
theorem weilPairing_self [CharZero F] [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n) (P : W.Point) (hP : n • P = 0) :
    weilPairing W n P P = 1 :=
  ABC3.Found.GaloisRep.weilPairingVal_alt W n hn P hP

/-- ★★★★**非退化性**——`P ≠ 0` なら `e_n(P, ·)` は自明でない。 -/
theorem weilPairing_nondeg (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (P : W.Point) (hP : n • P = 0) (hP0 : P ≠ 0) :
    ∃ Q : W.Point, n • Q = 0 ∧ weilPairing W n P Q ≠ 1 := by
  sorry

/-- ★★★★**Galois 同変性** `σ(e_n(P,Q)) = e_n(σP, σQ)`。

★★これが `det ρ = 円分指標`(`WeilPairing.lean` の `det_galRep_eq_cyclotomic`)を出す。 -/
theorem weilPairing_galois {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (n : ℕ) (σ : L ≃ₐ[K] L)
    (P Q : (W.baseChange L).toAffine.Point) :
    σ (weilPairing (W.baseChange L).toAffine n P Q)
      = weilPairing (W.baseChange L).toAffine n
          (ABC3.Interface.GaloisRep.galPoint W σ P)
          (ABC3.Interface.GaloisRep.galPoint W σ Q) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def weilPairing.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対 e_n の定義)",
    sectionId := "genell-thm-3-8" }

def weilPairing.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★★**2026-08-20: 定義は第 178 ブロック `weilPairingVal` になった**。`WeilSpec` という述語で witness の存在を言い、`Classical.choose` で値を取る。第 178 の `weilSpec_unique` が well-defined を保証する(0 ブロック)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——g_P の n 乗根の取り出し)" 19,
    .implicitStep
      "★`τ_Q(g_P)/g_P` が **`F` の元である**(定数関数である)ことは第 176 の `const_of_pow_eq_one` で済んだ。当初の見積もり 10-25 ブロックは第 176-178 の 3 ブロックで閉じた" 19,
    .implicitStep
      "★逸脱の記録: `[W.IsElliptic]` を仮定に足した(`weilPairingVal` の構成が群法則と生成点を使うため)。消費側 `det_galRep_eq_cyclotomic` は楕円曲線を扱うので後続に影響しない" 19 ]

def weilPairing_pow_eq_one.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対が 1 の n 乗根に値を取ること)",
    sectionId := "genell-thm-3-8" }

def weilPairing_pow_eq_one.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**2026-08-20: 第 179 ブロック `weilPairingVal_pow_eq_one` で証明された**。witness の `g^n = μ f_P` と `τ ∘ μ = μ`(第 168)から `(τg/g)^n = 1`(0 ブロック)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★`g_P(·+Q)^n = f_P([n]·+[n]Q) = f_P([n]·) = g_P(·)^n`——`[n]Q = 0` を使う 1 行だが、引き戻しの合成が可換であることが要る。第 168 の `aut_comp_mulByN` がそれである(0 ブロック)" 19,
    .implicitStep
      "★逸脱の記録: `[IsAlgClosed F]` を仮定に足した(第 137・第 139 が使う)。消費側は `[IsAlgClosed L]` の下なので後続に影響しない" 19 ]

def weilPairing_add_left.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_add_left.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★**2026-08-20: 第 184 + 第 191 ブロックで証明された**。3 点について**同じ `μ`・同じ `τ`** で witness を作り、第 183 の `(τg₁/g₁)(τg₂/g₂) = τg₃/g₃` を当てる。★`P₁ + P₂ = O` の場合は第 141 の `f_P·f_{−P} = c(x−x_P)^n` から別に出す(第 191)(0 ブロック)" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(a)(双線型性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★第 1 変数の加法性は `f_{P₁+P₂} = f_{P₁}·f_{P₂}·(主因子)` に帰着する。主因子の因子が 0 であることは第 182 の `elem_relation_of_add` で済んだ(0 ブロック)" 19,
    .implicitStep
      "★逸脱の記録: `[CharZero F]`・`[IsAlgClosed F]`・`[W.IsElliptic]` を仮定に足した。`μ` の存在(第 125)と `IsUnit (2:F)` に要る。消費側 `det_galRep_eq_cyclotomic` は `[CharZero K]`・`[IsAlgClosed L]` の下なので後続に影響しない" 19 ]

def weilPairing_self.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_self.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★★**2026-08-20: 第 190 + 第 191 ブロックで証明された**。`H := ∏_{i<n} τ_{iP'}^*(g)`(`nP' = P`)の `n` 乗が定数になり、`τ_{P'}^* H = H` の伸縮で `τ_P^*(g) = g`(0 ブロック)" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(b)(交代性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★`∏_{j} τ_{jP}(g_P)` の telescoping。★`n` 等分点 `P'` は第 186 で**分割多項式を使わずに**取れた(位数 `#E[n] = n²` の数え上げだけ)。★★極の位数 `−n` は第 187(値群の全射性)、積が定数になることは第 188、交換則 `τ_T ∘ [n]^* = [n]^* ∘ τ_{nT}` は第 189(0 ブロック)" 19,
    .implicitStep
      "★逸脱の記録: `hchar` の範囲が `k ≤ n²` に広がった(第 186 が `#E[n²] = n⁴` を使うため)。`[CharZero F]` の下では自動的に満たされる" 19 ]

def weilPairing_nondeg.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_nondeg.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(c)(非退化性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★★`e_n(P,·) ≡ 1` なら `g_P` が `[n]^*` の像に入り、`f_P` が `n` 乗になる——そこから `P = 0` を出す。Kummer 理論型の議論(20-50 ブロック)" 19,
    .implicitStep
      "★★★★2026-08-20 の測定: これは `F(E)/[n]^*F(E)` が **`E[n]` を Galois 群とする Galois 拡大**であることを要求する。★`τ_Q(g) = g`(∀Q ∈ E[n])⟹ `g ∈ [n]^*F(E)` の段である。★★mathlib の `IsGalois`・`FixedPoints` は使えるが、`E[n]` の作用が忠実かつ次数 `n²` であることを別途示す必要がある" 19,
    .implicitStep
      "★`E[n] ≃ (ℤ/n)²`(第 65-72 ブロック、**Found に済**)がここで消費される(0 ブロック)" 19 ]

def weilPairing_galois.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_galois.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(d)(Galois 同変性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★構成が関手的なら自動だが、明示式(`τ_Q(g_P)/g_P`)からだと `σ(f_P) = f_{σP}` と `σ(g_P) = g_{σP}` を別途示すことになる(10-25 ブロック)" 19,
    .implicitStep
      "★★★★2026-08-20 の見通し: 第 178 の `WeilSpec` は**データの存在**として書いてあるので、`σ` で witness を輸送すれば良い。★`σ : L ≃ₐ[K] L` は `W.CoordinateRing` と `W.FunctionField` にそれぞれ環同型を誘導する(係数が `K` にあるから)。★★`WeilSpec W n P Q c` の witness に `σ` を当てると `WeilSpec W n (σP) (σQ) (σc)` の witness になる——これが本節点の中身である" 19,
    .implicitStep
      "★`galPoint`(σ の点への作用)は **Interface に定義済**であり posit ではない(0 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep

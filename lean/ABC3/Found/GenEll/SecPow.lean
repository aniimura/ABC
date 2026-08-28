/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AmpleCover
import ABC3.Found.Arakelov.TrivTensor
import ABC3.Meta.Claim

/-!
# ★★★★★★★**`X_{s^{⊗k}} = X_s`** —— 段 E2 の後半（次数揃え）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★これは段 E2 の後半である

`§9-817` の有限被覆の切断は、それぞれ**別の** `M^{⊗n_j}` に属する。
★射 `X ⟶ ℙᴺ_R` を作るには次数を揃える必要があり、
`L = lcm(n_j)` として `s_j^{⊗(L/n_j)} ∈ Γ(X, M^{⊗L})` に置き換える。

★★そのとき**非消失軌跡が変わらない**ことが要る——それが本ファイルである:

    `X_{s^{⊗k}} = X_s`   （`k ≥ 1`）

## ★★★★★機構 —— 大域の `X_{s⊗t} = X_s ∩ X_t`

`TrivTensor.lean`（`§9-736` 頃）にあったのは**チャートを固定した形**
（`basicOpen_trivValue_tensor`）と片側の包含（`le_nonVanishing_tensor`）だけだった。
★本ファイルはそれを**大域の等式**に上げる（`nonVanishing_tensor`）。

★★機構は「各点で `A` と `B` が**同時に**自明化する開 `V` を取る」ことだけである
（`exists_common_triv`——2 つの被覆篩から `V_A ⊓ V_B` を取り、`trivialOfLe` で降ろす）。
★★★そこで `nonVanishing_inf`（段 D2）が両辺を `V` の中で `basicOpen` に直し、
`basicOpen_trivValue_tensor` が積を交わりに直す。

## ★★★冪への持ち上げ

`s^{⊗k}`（`secPow`）は `s ⊗ₜ s^{⊗(k-1)}` で定める（`presheafTensorPow` の再帰に合わせる）。
★`X_{s^{⊗1}} = X_s ⊓ X_{1} = X_s ⊓ ⊤ = X_s`（単位の切断 `1` は消えない）。
★★`k+1` の段は `X_s ⊓ X_s = X_s` で潰れる。

## ★残っている段（明示）

★★**`lcm` の帳簿は本ファイルに無い**——有限被覆の次数 `n_j` の最小公倍数を取り、
各 `j` で `k_j = L/n_j` として `s_j^{⊗k_j}` を作る段である。
★★★`X_{s^{⊗k}} = X_s` があるので被覆は保たれる（それが本ファイルの目的である）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★(1) 2 つが同時に自明化する開 -/

/-- ★★**各点で `A` と `B` が同時に自明化する開が取れる**。

★2 つの被覆篩から `V_A ⊓ V_B` を取り、`trivialOfLe` で降ろすだけである。 -/
theorem exists_common_triv {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) (x : X) :
    ∃ (V : X.Opens), x ∈ V ∧
      Nonempty ((restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V)) ∧
      Nonempty ((restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) := by
  obtain ⟨SA, hSA, htrivA⟩ := hA ⊤
  obtain ⟨VA, iA, hiA, hxA⟩ := hSA x trivial
  obtain ⟨eA⟩ := htrivA iA hiA
  obtain ⟨SB, hSB, htrivB⟩ := hB ⊤
  obtain ⟨VB, iB, hiB, hxB⟩ := hSB x trivial
  obtain ⟨eB⟩ := htrivB iB hiB
  exact ⟨VA ⊓ VB, ⟨hxA, hxB⟩,
    ⟨trivialOfLe A (inf_le_left : VA ⊓ VB ≤ VA) eA⟩,
    ⟨trivialOfLe B (inf_le_right : VA ⊓ VB ≤ VB) eB⟩⟩

/-! ## ★★★★★★★(2) 大域の `X_{s⊗t} = X_s ∩ X_t` -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**`X_{s⊗t} = X_s ∩ X_t`**（大域の形）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`TrivTensor.lean` はチャートを固定した形と片側の包含までだった。
★★本定理は「各点で同時に自明化する開を取る」ことでそれを大域の等式に上げる。 -/
theorem nonVanishing_tensor {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    nonVanishing (A ⊗ B) (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = nonVanishing A s ⊓ nonVanishing B t := by
  apply le_antisymm
  · intro x hx
    obtain ⟨V, hxV, ⟨eA⟩, ⟨eB⟩⟩ := exists_common_triv hA hB x
    have h1 := nonVanishing_inf (A ⊗ B) V (tensorTriv eA eB)
      (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
    have h2 := basicOpen_trivValue_tensor eA eB s t
    have h3 : x ∈ nonVanishing (A ⊗ B) (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t) ⊓ V := ⟨hx, hxV⟩
    rw [h1, h2] at h3
    have hA' := nonVanishing_inf A V eA s
    have hB' := nonVanishing_inf B V eB t
    rw [← hA', ← hB'] at h3
    exact ⟨h3.1.1, h3.2.1⟩
  · intro x hx
    obtain ⟨V, hxV, ⟨eA⟩, ⟨eB⟩⟩ := exists_common_triv hA hB x
    have hA' := nonVanishing_inf A V eA s
    have hB' := nonVanishing_inf B V eB t
    have h4 : x ∈ X.basicOpen (trivValue A V eA s) ⊓ X.basicOpen (trivValue B V eB t) := by
      rw [← hA', ← hB']
      exact ⟨⟨hx.1, hxV⟩, ⟨hx.2, hxV⟩⟩
    exact le_nonVanishing_tensor eA eB s t h4

/-! ## ★★★★★★★(3) 切断のテンソル冪 -/

/-- ★**切断のテンソル冪 `s^{⊗k}`**（`presheafTensorPow` の再帰に合わせる）。 -/
noncomputable def secPow (M : X.PresheafOfModules) (s : (M.obj (op ⊤) : Type)) :
    ∀ k : ℕ, ((presheafTensorPow M k).obj (op ⊤) : Type)
  | 0 => (1 : (Γ(X, (⊤ : X.Opens)) : Type))
  | k + 1 => s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] (secPow M s k)

/-- ★単位の切断 `1` はどこでも消えない。 -/
theorem nonVanishing_unit_one : nonVanishing (𝟙_ X.PresheafOfModules)
    ((1 : (Γ(X, (⊤ : X.Opens)) : Type))) = ⊤ := by
  rw [nonVanishing_unit]
  exact X.basicOpen_of_isUnit isUnit_one

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**`X_{s^{⊗k}} = X_s`**（`k ≥ 1`）——次数を揃えても被覆は変わらない。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが段 E2 の「次数揃え」を可能にする。 -/
theorem nonVanishing_secPow (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s : (M.obj (op ⊤) : Type)) : ∀ k : ℕ,
      nonVanishing (presheafTensorPow M (k + 1)) (secPow M s (k + 1)) = nonVanishing M s
  | 0 => by
      show nonVanishing (M ⊗ presheafTensorPow M 0)
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] secPow M s 0) = _
      rw [nonVanishing_tensor hM (isLocallyTrivial_presheafTensorPow hM 0)]
      show nonVanishing M s ⊓ nonVanishing (𝟙_ X.PresheafOfModules)
        ((1 : (Γ(X, (⊤ : X.Opens)) : Type))) = _
      rw [nonVanishing_unit_one, inf_top_eq]
  | k + 1 => by
      show nonVanishing (M ⊗ presheafTensorPow M (k + 1))
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] secPow M s (k + 1)) = _
      rw [nonVanishing_tensor hM (isLocallyTrivial_presheafTensorPow hM (k + 1)),
        nonVanishing_secPow M hM s k, inf_idem]

/-! ## ★出典の紐付け(`.src`) -/

def exists_common_triv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(各点で 2 つの層が同時に自明化する開が取れる)",
    sectionId := "genell-prop-1-4" }

def nonVanishing_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_{s⊗t} = X_s ∩ X_t——大域の形)",
    sectionId := "genell-prop-1-4" }

def secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——切断のテンソル冪 s^{⊗k})",
    sectionId := "genell-prop-1-4" }

def nonVanishing_secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_{s^{⊗k}} = X_s——次数を揃えても被覆は変わらない)",
    sectionId := "genell-prop-1-4" }

def nonVanishing_secPow.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "basicOpen_trivValue_tensor / le_nonVanishing_tensor(チャート固定の形)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.basicOpen_trivValue_tensor") 6,
    .citation "[ABC3]" "nonVanishing_inf(段 D2——X_s ⊓ V = basicOpen (trivValue))"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_inf") 6,
    .citation "[ABC3]" "isLocallyTrivial_presheafTensorPow(テンソル冪も局所自明)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isLocallyTrivial_presheafTensorPow") 6,
    .implicitStep
      ("★**lcm の帳簿は本ファイルに無い**——有限被覆の次数 n_j の最小公倍数 L を取り、" ++
       "各 j で k_j = L/n_j として s_j^{⊗k_j} を作る段である。" ++
       "★★X_{s^{⊗k}} = X_s があるので被覆は保たれる(それが本ファイルの目的である)") 6 ]

end ABC3.Found.GenEll

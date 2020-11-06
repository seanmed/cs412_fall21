from django import forms

class SearchSpotifyAPIForm(forms.Form):
    twitter_username= forms.CharField(label= 'twitter_username')
